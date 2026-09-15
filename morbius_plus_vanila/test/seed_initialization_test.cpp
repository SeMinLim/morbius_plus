#include "MorbiusPlus.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
using namespace std;


static void require( bool condition, const string &message ) {
	if ( condition ) return;
	cerr << "Seed initialization test failed: " << message << "\n";
	exit(1);
}

static Config testConfig( void ) {
	Config config = {};
	config.alphabetMode = ALPHABET_DNA;
	config.motifLength = 10;
	config.maxUpdateNum = 128;
	config.scoreThreshold = 1.0;
	config.randomSeed = 1;
	config.threadNum = 1;
	config.outputMotifNum = 16;
	return config;
}

// Exhaustive legal anchors are 0, 1, and 2 for this width-10 fixture.
// Count every column by its actual characters, independently of encoded BPMs.
static size_t optimalAnchor( const Config &config,
			     const Dataset &dataset,
			     const SeedModel &model ) {
	double bestRank = -numeric_limits<double>::infinity();
	size_t bestAnchor = 0;
	for ( size_t anchor = 0; anchor <= config.motifLength - model.seedLength; anchor ++ ) {
		vector<string> sites;
		for ( size_t sequenceIdx = 0; sequenceIdx < model.sampleNum; sequenceIdx ++ ) {
			const string &sequence = dataset.sequences[sequenceIdx];
			for ( size_t offset = 0; offset + config.motifLength <= sequence.size(); offset ++ ) {
				if ( sequence.compare(offset + anchor, model.seedLength, model.seedString) != 0 ) continue;
				sites.push_back(sequence.substr(offset, config.motifLength));
				break;
			}
		}
		if ( sites.empty() ) continue;
		size_t agreementSum = 0;
		for ( size_t column = 0; column < config.motifLength; column ++ ) {
			size_t maximum = 0;
			for ( char symbol : dataset.alphabet ) {
				size_t count = 0;
				for ( const string &site : sites ) if ( site[column] == symbol ) count ++;
				maximum = max(maximum, count);
			}
			agreementSum += maximum;
		}
		double rank = ((double)sites.size() / model.sampleNum) *
			      ((double)agreementSum / (sites.size() * config.motifLength) -
			       1.0 / dataset.alphabetSize);
		if ( rank > bestRank ) {
			bestRank = rank;
			bestAnchor = anchor;
		}
	}
	require(isfinite(bestRank), "fixture seed has no legal anchor");
	return bestAnchor;
}

static void checkSharedModels( const vector<SeedModel> &models ) {
	require(models.size() == 16, "missing shared seed assignments");
	const SeedModel &shared = models[0];
	for ( const SeedModel &model : models ) {
		require(model.valid == shared.valid && model.sampleNum == shared.sampleNum &&
			model.seedLength == shared.seedLength && model.seedCode == shared.seedCode &&
			model.seedString == shared.seedString && model.seedSupport == shared.seedSupport &&
			model.seedExpectedSupport == shared.seedExpectedSupport && model.seedRank == shared.seedRank &&
			model.anchorOffset == shared.anchorOffset && model.anchorRank == shared.anchorRank &&
			model.guidedOffsets == shared.guidedOffsets,
			"pipelines did not receive the same complete seed model");
	}
}

static void checkFallback( Config config, const string &sequence, size_t expectedSeedNum ) {
	Dataset dataset;
	configureAlphabet(config.alphabetMode, &dataset);
	dataset.sequenceLength = sequence.size();
	dataset.sequences.assign(4, sequence);
	dataset.names.assign(4, "fallback");
	config.maxUpdateNum = 1;
	validateWorkload(&config, &dataset);
	vector<SeedModel> models;
	buildSeedModels(&config, &dataset, models);
	checkSharedModels(models);
	size_t validNum = 0;
	for ( const SeedModel &model : models ) if ( model.valid ) validNum ++;
	require(validNum == expectedSeedNum, "shared seed validity or random fallback changed");
	vector<PipelineResult> results;
	runPipelines(&config, &dataset, models, results);
	SeedModel randomModel = {};
	for ( size_t pipelineIdx = 0; pipelineIdx < models.size(); pipelineIdx ++ ) {
		const vector<uint32_t> &offsets = results[pipelineIdx].initialOffsets;
		require(offsets.size() == dataset.sequences.size(), "fallback lost initial offsets");
		for ( uint32_t offset : offsets ) {
			require(offset + config.motifLength <= sequence.size(), "fallback offset out of bounds");
		}
		if ( models[pipelineIdx].valid == false ) {
			vector<uint32_t> expectedOffsets;
			initializeOffsets(&config, &dataset, &randomModel, (int)pipelineIdx, expectedOffsets);
			require(offsets == expectedOffsets, "invalid shared seed did not use existing random fallback");
		}
		if ( sequence.size() == config.motifLength ) {
			require(offsets == vector<uint32_t>(4, 0), "identical legal initial arrays were altered");
		}
	}
}

int main( int argc, char **argv ) {
	require(argc == 2, "DNA fixture argument is required");
	Config config = testConfig();
	Dataset dataset;
	configureAlphabet(config.alphabetMode, &dataset);
	readFASTA(argv[1], &dataset);
	validateWorkload(&config, &dataset);
	require(NUMPIPELINE == 16, "this experiment requires 16 pipelines");

	vector<SeedModel> models;
	buildSeedModels(&config, &dataset, models);
	checkSharedModels(models);

	// Golden shared seed and pipeline 0 offsets from the original single-seed
	// initialization in commit 7a04641d266f4ca8d921cccf965cae08b3b3192a.
	const vector<uint32_t> legacyInitialOffsets = {
		37, 46, 52, 48, 8, 49, 18, 4, 36, 41, 25, 8, 18, 10, 19, 11,
		46, 38, 33, 32, 30, 34, 5, 19, 0, 39, 34, 31, 44, 39, 15, 35
	};
	require(models[0].valid && models[0].seedCode == 23527, "original Markov seed winner changed");
	require(models[0].anchorOffset == optimalAnchor(config, dataset, models[0]),
		"shared anchor does not maximize the original anchor score");
	require(models[0].seedString == "CCGTTGCT" && models[0].seedSupport == 1 &&
		models[0].anchorOffset == 1, "legacy highest-ranked seed metadata changed");
	require(abs(models[0].seedExpectedSupport - 0.0029554144597797507) < 1.0e-14 &&
		abs(models[0].seedRank - 18.34110251631185) < 1.0e-12 &&
		models[0].anchorRank == 0.0234375, "legacy seed or anchor scoring changed");

	vector<PipelineResult> results;
	runPipelines(&config, &dataset, models, results);
	require(results[0].initialOffsets == legacyInitialOffsets, "pipeline 0 initialization changed");
	bool hasDifferentInitialOffsets = false;
	for ( size_t pipelineIdx = 0; pipelineIdx < models.size(); pipelineIdx ++ ) {
		const PipelineResult &result = results[pipelineIdx];
		require(result.bestOffsets.size() == dataset.sequences.size() &&
			result.bestStrands.size() == dataset.sequences.size(), "best state lost a sequence");
		for ( size_t seqIdx = 0; seqIdx < dataset.sequences.size(); seqIdx ++ ) {
			require(result.bestOffsets[seqIdx] + config.motifLength <= dataset.sequenceLength &&
				result.bestStrands[seqIdx] <= STRAND_REVERSE, "best site is out of bounds");
		}
		vector<uint32_t> bestBPM;
		buildBPM(&config, &dataset, result.bestOffsets, result.bestStrands, bestBPM);
		require(calculateAgreementScore(&config, &dataset, bestBPM) == result.bestScore,
			"best score does not match stored offsets and strands");
		vector<uint32_t> expectedOffsets;
		initializeOffsets(&config, &dataset, &models[0], (int)pipelineIdx, expectedOffsets);
		require(result.initialOffsets == expectedOffsets,
			"shared seed initialization did not preserve the pipeline's own random stream");
		if ( result.initialOffsets != results[0].initialOffsets ) hasDifferentInitialOffsets = true;
		for ( uint32_t offset : expectedOffsets ) {
			require(offset + config.motifLength <= dataset.sequenceLength, "initial offset out of bounds");
		}
	}
	require(hasDifferentInitialOffsets, "fixture cannot detect accidental copying of pipeline 0 offsets");

	// Thread scheduling and requested output count must not affect discovery.
	config.threadNum = 4;
	config.outputMotifNum = 1;
	vector<PipelineResult> repeatedResults;
	runPipelines(&config, &dataset, models, repeatedResults);
	for ( size_t pipelineIdx = 0; pipelineIdx < results.size(); pipelineIdx ++ ) {
		const PipelineResult &first = results[pipelineIdx];
		const PipelineResult &second = repeatedResults[pipelineIdx];
		require(first.initialOffsets == second.initialOffsets && first.bestOffsets == second.bestOffsets &&
			first.bestStrands == second.bestStrands &&
			first.bestScore == second.bestScore && first.updateNum == second.updateNum &&
			first.thresholdReached == second.thresholdReached,
			"thread count or output motif count changed a pipeline result");
	}

	config.motifLength = 8;
	checkFallback(config, "ACGTACGT", 16);
	config.motifLength = 10;
	checkFallback(config, "AAAAAAAAAAAA", 0);
	cout << "All shared-seed initialization tests passed.\n";
	return 0;
}
