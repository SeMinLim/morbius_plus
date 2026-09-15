#include "MorbiusPlus.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
using namespace std;

// Exercise the internal sampler and update loop without changing their visibility.
CandidateSite sampleCandidate( const Config *, const Dataset *, const string &,
			       const vector<double> &, RandomGenerator * );
double calculateCandidateLogProb( const Config *, const Dataset *, const string &,
				 size_t, uint8_t, const vector<double> & );
void runPipeline( const Config *, const Dataset *, const SeedModel *, int, uint64_t,
		  PipelineResult * );

static void require( bool condition, const string &message ) {
	if ( condition ) return;
	cerr << "Strand test failed: " << message << "\n";
	exit(1);
}

static Config testConfig( size_t width ) {
	Config config = {};
	config.alphabetMode = ALPHABET_DNA;
	config.motifLength = width;
	config.maxUpdateNum = 256;
	config.scoreThreshold = 1.0;
	config.randomSeed = 1;
	config.threadNum = 1;
	config.outputMotifNum = 16;
	return config;
}

static Dataset testDataset( const vector<string> &sequences ) {
	Dataset dataset;
	configureAlphabet(ALPHABET_DNA, &dataset);
	dataset.sequences = sequences;
	dataset.sequenceLength = sequences[0].size();
	dataset.names.assign(sequences.size(), "test");
	return dataset;
}

// Independent string construction keeps reconstruction separate from getSiteSymbol.
static string siteString( const string &sequence, uint32_t offset, uint8_t strand,
			  size_t width ) {
	string site = sequence.substr(offset, width);
	if ( strand == STRAND_REVERSE ) {
		reverse(site.begin(), site.end());
		for ( char &base : site ) {
			switch ( base ) {
				case 'A': base = 'T'; break;
				case 'C': base = 'G'; break;
				case 'G': base = 'C'; break;
				case 'T': base = 'A'; break;
				default: require(false, "invalid test base");
			}
		}
	}
	return site;
}

static vector<uint32_t> reconstructBPM( const Config &config, const Dataset &dataset,
				     const vector<uint32_t> &offsets,
				     const vector<uint8_t> &strands ) {
	vector<uint32_t> bpm(4 * config.motifLength, PSEUDOCOUNT);
	for ( size_t seqIdx = 0; seqIdx < dataset.sequences.size(); seqIdx ++ ) {
		string site = siteString(dataset.sequences[seqIdx], offsets[seqIdx], strands[seqIdx],
					config.motifLength);
		for ( size_t column = 0; column < config.motifLength; column ++ ) {
			bpm[dataset.alphabet.find(site[column]) * config.motifLength + column] ++;
		}
	}
	return bpm;
}

static uint64_t agreement( const vector<uint32_t> &bpm, size_t width ) {
	uint64_t score = 0;
	for ( size_t column = 0; column < width; column ++ ) {
		uint32_t maximum = 0;
		for ( size_t symbol = 0; symbol < 4; symbol ++ ) maximum = max(maximum, bpm[symbol * width + column]);
		score += maximum;
	}
	return score;
}

static void checkJointSampling( void ) {
	Config config = testConfig(1);
	Dataset dataset = testDataset({"AC"});
	vector<double> lpm = {0.0, 1.0, 2.0, 3.0};
	RandomGenerator random;
	initializeRandomGenerator(&random, 1947);
	const size_t sampleNum = 60000;
	vector<size_t> counts(4, 0);
	for ( size_t trial = 0; trial < sampleNum; trial ++ ) {
		CandidateSite site = sampleCandidate(&config, &dataset, "AC", lpm, &random);
		require(site.offset < 2 && site.strand <= STRAND_REVERSE, "joint sample out of bounds");
		counts[(size_t)site.strand * 2 + site.offset] ++;
	}
	// Integer log2 values give exact PWL weights. Both orientation totals belong
	// to the same normalization: forward 3/15, reverse 12/15, not one half each.
	const double expected[] = {1.0 / 15, 2.0 / 15, 8.0 / 15, 4.0 / 15};
	for ( size_t idx = 0; idx < counts.size(); idx ++ ) {
		require(abs((double)counts[idx] / sampleNum - expected[idx]) < 0.01,
			"offset/strand probabilities are not jointly normalized");
	}

	// Put the orientations in separate CDF segments, with unequal total mass.
	string sequence(32, 'A');
	dataset = testDataset({sequence});
	size_t reverseNum = 0;
	for ( size_t trial = 0; trial < sampleNum; trial ++ ) {
		CandidateSite site = sampleCandidate(&config, &dataset, sequence, lpm, &random);
		require(site.offset < sequence.size() && site.strand <= STRAND_REVERSE,
			"segmented sample out of bounds");
		if ( site.strand == STRAND_REVERSE ) reverseNum ++;
	}
	require(abs((double)reverseNum / sampleNum - 8.0 / 9) < 0.01,
		"separate segments lost the joint orientation mass");
}

static void checkOrientedCounts( void ) {
	Config config = testConfig(4);
	Dataset dataset = testDataset({"ACGAAG", "GGTCGT"});
	vector<uint32_t> offsets = {0, 2};
	vector<uint8_t> strands = {STRAND_FORWARD, STRAND_REVERSE};
	vector<uint32_t> expected = reconstructBPM(config, dataset, offsets, strands);
	vector<uint32_t> bpm;
	buildBPM(&config, &dataset, offsets, strands, bpm);
	require(bpm == expected, "BPM does not align reverse-complement sites");
	vector<uint32_t> count;
	buildResultCount(&config, &dataset, offsets, strands, count);
	for ( size_t idx = 0; idx < count.size(); idx ++ ) {
		require(count[idx] + PSEUDOCOUNT == expected[idx], "output count differs from oriented BPM");
	}
	require(buildConsensus(&config, &dataset, count) == "ACGA", "known oriented consensus is wrong");
	for ( size_t column = 0; column < config.motifLength; column ++ ) {
		uint32_t sum = 0;
		for ( size_t symbol = 0; symbol < 4; symbol ++ ) sum += count[symbol * config.motifLength + column];
		require(sum == dataset.sequences.size(), "a sequence contributes twice to the PWM");
	}
	vector<double> lpm(count.size(), 0.0);
	for ( size_t idx = 0; idx < lpm.size(); idx ++ ) lpm[idx] = (double)idx;
	double expectedLogProb = lpm[0] + lpm[5] + lpm[10] + lpm[3];
	require(calculateCandidateLogProb(&config, &dataset, "GGTCGT", 2, STRAND_REVERSE, lpm) ==
		expectedLogProb, "candidate evaluation uses the wrong reversed column or complement");
}

static void checkIncrementalUpdates( const vector<string> &sequences, size_t width ) {
	Config config = testConfig(width);
	Dataset dataset = testDataset(sequences);
	SeedModel randomModel = {};
	vector<uint32_t> offsets;
	initializeOffsets(&config, &dataset, &randomModel, 0, offsets);
	vector<uint8_t> strands(sequences.size(), STRAND_FORWARD);
	vector<uint32_t> expectedOffsets = offsets;
	vector<uint8_t> expectedStrands = strands;
	uint64_t bestScore = agreement(reconstructBPM(config, dataset, offsets, strands), width);
	RandomGenerator random;
	initializeRandomGenerator(&random, config.randomSeed ^ 0xd2b74407b1ce6e93ULL);
	size_t reverseSamples = 0;
	for ( uint64_t updateIdx = 0; updateIdx < config.maxUpdateNum; updateIdx ++ ) {
		size_t seqIdx = updateIdx % sequences.size();
		vector<uint32_t> bpm = reconstructBPM(config, dataset, offsets, strands);
		string oldSite = siteString(sequences[seqIdx], offsets[seqIdx], strands[seqIdx], width);
		for ( size_t column = 0; column < width; column ++ ) {
			bpm[dataset.alphabet.find(oldSite[column]) * width + column] --;
		}
		vector<double> lpm;
		buildLPM(bpm, lpm);
		CandidateSite site = sampleCandidate(&config, &dataset, sequences[seqIdx], lpm, &random);
		offsets[seqIdx] = site.offset;
		strands[seqIdx] = site.strand;
		if ( site.strand == STRAND_REVERSE ) reverseSamples ++;
		uint64_t score = agreement(reconstructBPM(config, dataset, offsets, strands), width);
		if ( score > bestScore ) {
			bestScore = score;
			expectedOffsets = offsets;
			expectedStrands = strands;
		}
	}
	require(reverseSamples > 0, "replay fixture did not sample reverse sites");
	if ( sequences.size() > 1 ) {
		require(find(expectedStrands.begin(), expectedStrands.end(), STRAND_REVERSE) != expectedStrands.end(),
			"replay fixture has no reverse site in its best state");
		require(expectedStrands != strands, "replay fixture must distinguish best from final strands");
	}
	PipelineResult result;
	// An unreachable threshold forces removal/addition across every update,
	// including when the only sequence is also the next sequence to remove.
	runPipeline(&config, &dataset, &randomModel, 0, numeric_limits<uint64_t>::max(), &result);
	require(result.bestOffsets == expectedOffsets && result.bestStrands == expectedStrands &&
		result.bestScore == bestScore, "incremental BPM/LPM or saved best strands differ from full reconstruction");
	require(result.updateNum == config.maxUpdateNum && result.thresholdReached == false,
		"update-limit termination is incorrect");
	runPipeline(&config, &dataset, &randomModel, 0, 0, &result);
	require(result.updateNum == 0 && result.thresholdReached &&
		result.bestOffsets == result.initialOffsets &&
		result.bestStrands == vector<uint8_t>(sequences.size(), STRAND_FORWARD),
		"initial threshold stop altered the forward initial state");
}

int main( void ) {
	checkJointSampling();
	checkOrientedCounts();
	checkIncrementalUpdates({"AACGATT", "CCATGCT", "GCTTTAG", "TGGCAAC", "TTAGCAA"}, 5);
	checkIncrementalUpdates({"ACGCAA"}, 6);
	cout << "All strand sampling and state tests passed.\n";
	return 0;
}
