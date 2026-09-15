#include "Refinement.h"

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
	cerr << "Refinement test failed: " << message << "\n";
	exit(1);
}

static void closeTo( double actual, double expected, const string &message ) {
	require(abs(actual - expected) < 1e-10, message);
}

static Config testConfig( size_t width ) {
	Config config = {};
	config.alphabetMode = ALPHABET_DNA;
	config.motifLength = width;
	config.outputMotifNum = NUMPIPELINE;
	return config;
}

static Dataset testDataset( const vector<string> &sequences ) {
	Dataset dataset;
	configureAlphabet(ALPHABET_DNA, &dataset);
	dataset.sequences = sequences;
	dataset.sequenceLength = sequences[0].size();
	for ( size_t idx = 0; idx < sequences.size(); idx ++ ) dataset.names.push_back("sequence_" + to_string(idx));
	return dataset;
}

static string reverseComplement( string sequence ) {
	reverse(sequence.begin(), sequence.end());
	const string alphabet = "ACGT";
	for ( char &base : sequence ) base = alphabet[3 - alphabet.find(base)];
	return sequence;
}

static string siteString( const string &sequence, size_t offset, uint8_t strand, size_t width ) {
	string site = sequence.substr(offset, width);
	return strand == STRAND_REVERSE ? reverseComplement(site) : site;
}

// Small integer combinations give an independent oracle for every 2x2 table.
static double choose( size_t n, size_t k ) {
	if ( k > n ) return 0.0;
	double result = 1.0;
	for ( size_t idx = 1; idx <= k; idx ++ ) result *= (double)(n + 1 - idx) / idx;
	return result;
}

static double exactFisher( size_t primaryHits, size_t controlHits, size_t primaryNum, size_t controlNum ) {
	size_t hits = primaryHits + controlHits;
	double sum = 0.0;
	for ( size_t idx = primaryHits; idx <= min(primaryNum, hits); idx ++ ) {
		sum += choose(primaryNum, idx) * choose(controlNum, hits - idx);
	}
	return sum / choose(primaryNum + controlNum, hits);
}

static void checkFisherAndTies( void ) {
	for ( size_t primaryNum = 1; primaryNum <= 7; primaryNum ++ ) {
		for ( size_t controlNum = 1; controlNum <= 7; controlNum ++ ) {
			for ( size_t primaryHits = 0; primaryHits <= primaryNum; primaryHits ++ ) {
				for ( size_t controlHits = 0; controlHits <= controlNum; controlHits ++ ) {
					closeTo(exp(fisherLogPvalue(primaryHits, controlHits, primaryNum, controlNum)),
						exactFisher(primaryHits, controlHits, primaryNum, controlNum), "one-sided Fisher table");
				}
			}
		}
	}
	vector<RefinementSite> primary = {{6, 0, 0}, {5, 0, 0}, {5, 0, 0}, {5, 0, 0}};
	vector<RefinementSite> control = {{5, 0, 0}, {5, 0, 0}, {1, 0, 0}, {1, 0, 0}};
	RefinementSelection selection = selectRefinementThreshold(primary, control);
	require(selection.valid && selection.scoreThreshold == 5.0, "best threshold at a tied score");
	require(selection.primarySiteNum == 4 && selection.controlSiteNum == 2,
		"threshold must include the complete Primary and Control score tie");
	closeTo(exp(selection.logPvalue), 15.0 / 70.0, "known tied-cutoff Fisher value");
	reverse(primary.begin(), primary.end());
	reverse(control.begin(), control.end());
	RefinementSelection reversed = selectRefinementThreshold(primary, control);
	closeTo(reversed.logPvalue, selection.logPvalue, "input ordering must not break score ties");
	require(!selectRefinementThreshold(primary, primary).valid,
		"equal Primary/Control hit fractions provide no enrichment");
	require(!selectRefinementThreshold({{0, 0, 0}, {-1, 0, 0}}, {{-2, 0, 0}, {-3, 0, 0}}).valid,
		"nonpositive LLR sites cannot support refinement");
	double extreme = fisherLogPvalue(10000, 0, 10000, 10000);
	require(isfinite(extreme) && extreme < -10000.0,
		"large enrichment must remain finite in log space even when p underflows");
}

static size_t code( const string &tuple ) {
	size_t result = 0;
	for ( char base : tuple ) result = 4 * result + string("ACGT").find(base);
	return result;
}

// Count concrete strings rather than reproducing the rolling-code implementation.
static vector<vector<double>> manualBackground( const Dataset &control ) {
	vector<vector<double>> probabilities;
	for ( size_t length = 1; length <= 3; length ++ ) {
		size_t tupleNum = (size_t)pow(4.0, (double)length);
		vector<double> counts(tupleNum, 1.0 / tupleNum);
		for ( const string &sequence : control.sequences ) {
			for ( const string &oriented : {sequence, reverseComplement(sequence)} ) {
				for ( size_t offset = 0; offset + length <= oriented.size(); offset ++ ) {
					counts[code(oriented.substr(offset, length))] += 1.0;
				}
			}
		}
		vector<double> conditional(tupleNum);
		for ( size_t start = 0; start < tupleNum; start += 4 ) {
			double sum = counts[start] + counts[start + 1] + counts[start + 2] + counts[start + 3];
			for ( size_t base = 0; base < 4; base ++ ) conditional[start + base] = counts[start + base] / sum;
		}
		probabilities.push_back(conditional);
	}
	return probabilities;
}

static double manualBackgroundScore( const string &site, const vector<vector<double>> &background ) {
	double score = 0.0;
	for ( size_t column = 0; column < site.size(); column ++ ) {
		size_t order = min((size_t)2, column);
		score += log2(background[order][code(site.substr(column - order, order + 1))]);
	}
	return score;
}

static RefinementSite manualBestSite( const string &sequence, size_t width,
		const vector<double> &pwm, const vector<vector<double>> &background ) {
	RefinementSite best = {-numeric_limits<double>::infinity(), 0, STRAND_FORWARD};
	for ( uint8_t strand = STRAND_FORWARD; strand <= STRAND_REVERSE; strand ++ ) {
		for ( size_t offset = 0; offset + width <= sequence.size(); offset ++ ) {
			string site = siteString(sequence, offset, strand, width);
			double score = -manualBackgroundScore(site, background);
			for ( size_t column = 0; column < width; column ++ ) {
				score += log2(pwm[string("ACGT").find(site[column]) * width + column]);
			}
			if ( score > best.score ) best = {score, (uint32_t)offset, strand};
		}
	}
	return best;
}

static void checkBackgroundAndStrands( void ) {
	Config config = testConfig(4);
	Dataset control = testDataset({"AAAAACGT", "CCAAGTTA", "CGCGGGTA"});
	RefinementBackground background;
	buildRefinementBackground(&control, &background);
	vector<vector<double>> expected = manualBackground(control);
	require(background.logProbability.size() == 3, "order-two background dimensions");
	for ( size_t order = 0; order < expected.size(); order ++ ) {
		require(background.logProbability[order].size() == expected[order].size(), "conditional table dimensions");
		for ( size_t idx = 0; idx < expected[order].size(); idx ++ ) {
			closeTo(exp2(background.logProbability[order][idx]), expected[order][idx],
				"Control-only, reverse-complement-symmetrized tuple probability");
		}
	}
	for ( size_t base = 0; base < 4; base ++ ) closeTo(background.letterFrequency[base], expected[0][base], "Control letter frequency");
	Dataset dataset = testDataset({"TTACGATT", "GGGTCGTT"});
	vector<double> pwm(16, 0.02);
	for ( size_t column = 0; column < 4; column ++ ) pwm[string("ACGT").find(string("ACGA")[column]) * 4 + column] = 0.94;
	vector<RefinementSite> scanned;
	scanRefinementSites(&config, &dataset, &background, pwm, scanned);
	require(scanned.size() == 2, "one best site per original sequence");
	for ( size_t seqIdx = 0; seqIdx < dataset.sequences.size(); seqIdx ++ ) {
		const string &sequence = dataset.sequences[seqIdx];
		for ( uint8_t strand = 0; strand < 2; strand ++ ) {
			for ( uint32_t offset = 0; offset + 4 <= sequence.size(); offset ++ ) {
				closeTo(calculateRefinementBackgroundLogProbability(&config, &dataset, &background, sequence, offset, strand),
					manualBackgroundScore(siteString(sequence, offset, strand, 4), expected),
					"background context must restart inside each oriented window");
			}
		}
		RefinementSite best = manualBestSite(sequence, 4, pwm, expected);
		closeTo(scanned[seqIdx].score, best.score, "manual likelihood ratio");
		require(scanned[seqIdx].offset == best.offset && scanned[seqIdx].strand == best.strand,
			"best-site orientation and original-coordinate offset");
	}
	require(scanned[0].offset == 2 && scanned[0].strand == STRAND_FORWARD &&
		scanned[1].offset == 3 && scanned[1].strand == STRAND_REVERSE, "planted forward and reverse sites");
	RefinementBackground uniform;
	for ( size_t order = 0; order < 3; order ++ ) uniform.logProbability.push_back(vector<double>((size_t)pow(4.0, order + 1), -2.0));
	scanRefinementSites(&config, &dataset, &uniform, vector<double>(16, 0.25), scanned);
	for ( const RefinementSite &site : scanned ) require(site.offset == 0 && site.strand == STRAND_FORWARD && site.score == 0.0,
		"equal site scores must choose forward, lowest offset");
}

static void checkSubsetRefinement( void ) {
	Config config = testConfig(8);
	// A fixed mixed-strand fixture whose second accepted scan re-admits four
	// sequences omitted by the first threshold. No Gibbs run is needed here.
	vector<string> primarySequences = {
		"ATAATACGATCGAAGCGCAAACTCGTAGCCAA", "CCGGATTCGATCGTACCCATAAGGGTCAGGAA",
		"CCTCTTCACGATCGATAGTCTGAGAGCGCGTA", "GATAAGGGTCGATCGTGCGTCCCCCGGTTTAG",
		"GTCTCGTCAACGATCGAGTCGTTCTCTTTTAA", "TTAGCGAGGATCGATCGTGAGTAAACGTAAAT",
		"TCCTTCTATCGACGATCGACTTCGCGGTGCTG", "GATTGATTATAGTCGATCGTCACTCATCGAAG",
		"TCAGGAGGGATTAACGATCGACGGGGTCCGGA", "TACCCATCACGGGTTCGATCGTTACTATCAAA",
		"GTCGGAGTAAAGACGACGATCGAACTGCTAAG", "AGCCATCGATCGTTCGCCGTTGTAGGCCCGGG",
		"TCCAAGACGATCGAGGGGCTAACCAGACGAGT", "CCGGGTGTCGATCGTTTCATCACAATTTCATG",
		"ACGCCGAAACGATCGACAAAATAAATGAACCC", "CAGCGAAAATCGATCGTCGCCATGTATTCAAT",
		"GGTTTTAAGCGATATTATGCTGATGTACAATA", "ACGCCCATCGTTGCCTCTTCTCTCGAGATTAT",
		"TTAATAACCCAGTCGCAAGGGACCTAATTTTC", "CCGGTGTTAACCATGATTGCCCTTGAGAGGCC",
		"AACTGCTATTACGCCGACTTATGGCGGTATTA", "GACAGATGCATACAAGTATATGGAACGTACGC",
		"GAAAAAAGAACTATAAGCGGATTAGGCAACGG", "ACACCCACCGAAGGCACCACCACCCAGACCTC",
		"ATGACACCACAAACCATCAGCCCACCCAGACT", "TACTGAACATATTACCGCAGTCGGTGCGATTA",
		"CGTGATATCTTTCTCACAATGGACGACCTCGA", "TGAATCCAGTCAAGCAGTGGAGATCAGCGGTA",
		"AGGAGAAAAAACGGTTGCTGCCCACCTACCGA", "GCGTGTTGACTGACATCAGGAACGTATCTAGA",
		"CGCGTGTTTCTGATTAACTCTCTTAGTCAATT", "CATGATAGCCTCACTGATACCGGACGGTGGAG",
		"TTTGATTACCATGAAACCATCGCATGAAGGCA", "TCGCGCACGTCCACCGGGATTTGGTACCCCCT",
		"TTCTAGGCAGGCACATCGCCAGGTGATTAGTT", "TATGAATCTATTGTCTTCTACGGTATTTATGA"
	};
	vector<string> controlSequences = {
		"TCTTGAGTCGGTAACCTCTCTCTCGCTCTATC", "GCATAGCCGCACTAAGTGATCCCTAAAAGGGA",
		"GAACGCGGGTAATACACAGCGGTCACCCGCGG", "AACTTAAAAAAAATAGACCCGCTCCGCATAGC",
		"AATAGAACGTCACTAGCCGACACGTCAACTAA", "ATCCTCGAATCGAGCCACATGGGAGGTGGGTG",
		"CCGTTGATGGATTAATTCTCAGTCGTGAAGTT", "CAATTTAGTAAGTCATATCCATGAAAGCCTTA",
		"AGGGACTCCCCTCTAAAGGTTTTCAATCCACT", "GATGACTTATGGGGTAGCGAGAAAAAAGTCGT",
		"CCGTATAGCTGACTGCGTCAAGTTGCTCAGCA", "CCGGTCAGGACGTACAAGTTCAATGAATCTGA",
		"GGGCCATGATCCGTAGGGCGCTAACTGGACAC", "GTTCGTCCTGCTGTACGTAGAGCCTCCGGACG",
		"AGGCACGCGGCGCGCACTAACGGCGTCTTTGC", "GCAGCCAGCTATCACTCGCTAACCATGCGATA",
		"GGAGCGCGCGACAGGAGTAGTGCGTACCTTCG", "ACAAAACGATTGCGATTCAACAGTAAGGCACT",
		"GGTAGCGAACTGGTACTACGGCCACCTCCAGG", "AGCTTGGGCTGTTATAGACATGGACCCATGTA",
		"CATGTTCGTTGTGATACGTTAGAGGGTACGAT", "TGTATTCACCCTTGTAGGTGAGAGAGTGGTCA",
		"CCCTTGCAAAGTCCGAACTGTGGTTTAACTAG", "GGTGTACGCCGCTACAAGCATCTATTCCTAGA",
		"GATCTGCTGCTACTACATGCTTCTATCTGGGG", "GGTACGGTTCTAGCACGCTTGCTAAAGGGCTT",
		"AGCCTAAATGCTATCGGGACCTCTCACTTGTG", "CCTGTCTGTGCCAATGCAATATGGCGATTATG",
		"TTTCAGCTCGTCAGTCGACTCGCCTTCACAAT", "TCTACTACACAAGTGGAGAAGCGAGATTGTAT",
		"CCCTATACATCCCTTGACGAAGGGCACCCACT", "AATCCTTCGATAATTTATCGGTTAGGAGGTGC",
		"GACCATCCTCAGAGCGGAGTTGCATGCCCAGT", "CAGGGTTCCCTGGGGGTCAATAAAGTATTTCG",
		"AGACGACGAAGGAAGACGGAAGCCGCTAAATG", "CGGAGGCTAGCAAGATGGAAACCGTCTCATGT"
	};
	Dataset primary = testDataset(primarySequences), control = testDataset(controlSequences);
	RefinementBackground background;
	buildRefinementBackground(&control, &background);
	vector<vector<double>> expectedBackground = manualBackground(control);
	OutputMotif initial{};
	initial.pipelineIdx = 12;
	initial.siteNum = primary.sequences.size();
	initial.sitePresent.assign(initial.siteNum, 1);
	initial.offsets = {5,6,7,8,9,10,11,12,13,14,15,5,6,7,8,9,8,3,8,12,13,20,5,3,8,20,15,5,24,17,9,23,15,9,23,6};
	initial.strands = {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,1,1};
	initial.count.assign(32, 0);
	for ( size_t idx = 0; idx < initial.siteNum; idx ++ ) {
		string site = siteString(primary.sequences[idx], initial.offsets[idx], initial.strands[idx], 8);
		for ( size_t column = 0; column < 8; column ++ ) initial.count[string("ACGT").find(site[column]) * 8 + column] ++;
	}
	RefinementResult result;
	refineMotif(&config, &primary, &control, &background, initial, &result);
	require(result.selection.valid && result.motif.refined && result.motif.pipelineIdx == 12, "candidate identity and refinement status");
	require(result.motif.siteNum == 20 && result.fittingSelection.primarySiteNum == 20,
		"accepted fitting sites must exclude unsupported Primary sequences");
	require(result.iterations > 0 && result.iterations <= REFINEMENTMAXITERATIONS && result.acceptedIterations > 0,
		"refinement stops within its separate iteration limit");
	require(result.acceptedLogPvalues.size() == result.acceptedIterations, "accepted-objective history length");
	require(result.acceptedIterations == 2 && result.acceptedLogPvalues[0] < result.initialLogPvalue,
		"first accepted PWM must improve the evaluated starting PWM");
	for ( size_t idx = 1; idx < result.acceptedLogPvalues.size(); idx ++ ) {
		require(result.acceptedLogPvalues[idx] < result.acceptedLogPvalues[idx - 1], "accepted enrichment must strictly improve");
	}
	vector<uint32_t> counts(32, 0);
	size_t primaryHits = 0, controlHits = 0;
	for ( size_t seqIdx = 0; seqIdx < primary.sequences.size(); seqIdx ++ ) {
		RefinementSite site = manualBestSite(primary.sequences[seqIdx], 8, result.scoringPWM, expectedBackground);
		bool present = site.score >= result.fittingSelection.scoreThreshold;
		require(result.motif.sitePresent[seqIdx] == present, "site mask must match the saved scoring PWM and threshold");
		if ( !present ) continue;
		primaryHits ++;
		require(result.motif.offsets[seqIdx] == site.offset && result.motif.strands[seqIdx] == site.strand,
			"reported positions must come from the accepted refinement scan");
		string selected = siteString(primary.sequences[seqIdx], site.offset, site.strand, 8);
		for ( size_t column = 0; column < 8; column ++ ) counts[string("ACGT").find(selected[column]) * 8 + column] ++;
	}
	for ( const string &sequence : control.sequences ) {
		if ( manualBestSite(sequence, 8, result.scoringPWM, expectedBackground).score >= result.fittingSelection.scoreThreshold ) controlHits ++;
	}
	require(counts == result.motif.count && primaryHits == result.motif.siteNum && controlHits == result.fittingSelection.controlSiteNum,
		"only oriented passing Primary sites contribute counts");
	closeTo(exp(result.fittingSelection.logPvalue), exactFisher(primaryHits, controlHits, 36, 36), "saved fitting selection enrichment");
	vector<double> finalPWM(32), initialPWM(32);
	for ( size_t idx = 0; idx < 32; idx ++ ) {
		finalPWM[idx] = (counts[idx] + 1.0) / (result.motif.siteNum + 4.0);
		initialPWM[idx] = (initial.count[idx] + 1.0) / (initial.siteNum + 4.0);
	}
	vector<RefinementSite> initialPrimary, initialControl;
	for ( const string &sequence : primary.sequences ) initialPrimary.push_back(manualBestSite(sequence, 8, initialPWM, expectedBackground));
	for ( const string &sequence : control.sequences ) initialControl.push_back(manualBestSite(sequence, 8, initialPWM, expectedBackground));
	RefinementSelection initialSelection = selectRefinementThreshold(initialPrimary, initialControl);
	size_t reentered = 0;
	for ( size_t idx = 0; idx < initialPrimary.size(); idx ++ ) {
		if ( initialPrimary[idx].score < initialSelection.scoreThreshold && result.motif.sitePresent[idx] ) reentered ++;
	}
	require(reentered == 4, "previously excluded sequences must be rescanned and allowed to re-enter");
	primaryHits = controlHits = 0;
	for ( const string &sequence : primary.sequences ) {
		if ( manualBestSite(sequence, 8, finalPWM, expectedBackground).score >= result.selection.scoreThreshold ) primaryHits ++;
	}
	for ( const string &sequence : control.sequences ) {
		if ( manualBestSite(sequence, 8, finalPWM, expectedBackground).score >= result.selection.scoreThreshold ) controlHits ++;
	}
	require(primaryHits == result.selection.primarySiteNum && controlHits == result.selection.controlSiteNum,
		"final enrichment counts must come from rescanning the final output PWM");
	closeTo(exp(result.selection.logPvalue), exactFisher(primaryHits, controlHits, 36, 36), "evaluated final PWM enrichment");
	for ( size_t column = 0; column < 8; column ++ ) {
		double probabilitySum = 0.0;
		uint32_t columnCount = 0;
		for ( size_t base = 0; base < 4; base ++ ) {
			columnCount += counts[base * 8 + column];
			probabilitySum += (counts[base * 8 + column] + 1.0) / (result.motif.siteNum + 4.0);
		}
		require(columnCount == 20, "excluded sequences must not contribute a site");
		closeTo(probabilitySum, 1.0, "pseudocount normalization uses participating site count");
	}
	RefinementResult unsupported;
	refineMotif(&config, &primary, &primary, &background, initial, &unsupported);
	require(!unsupported.selection.valid && unsupported.motif.siteNum == 0 && unsupported.acceptedIterations == 0,
		"identical Primary and Control must not fabricate a supported PWM");
}

static void checkRefinedDeduplication( void ) {
	Config config = testConfig(4);
	vector<PipelineResult> pipelines(4);
	vector<RefinementResult> results(4);
	const size_t siteNums[] = {4, 12, 4, 4};
	const uint32_t baseCounts[][4] = {{2, 1, 1, 0}, {5, 3, 3, 1}, {1, 2, 1, 0}, {4, 0, 0, 0}};
	for ( size_t idx = 0; idx < 4; idx ++ ) {
		pipelines[idx].bestScore = idx == 1 ? 100 : (10 - idx);
		results[idx].selection.valid = idx != 3;
		results[idx].motif.pipelineIdx = (int)idx;
		results[idx].motif.refined = true;
		results[idx].motif.siteNum = siteNums[idx];
		for ( size_t base = 0; base < 4; base ++ ) {
			for ( size_t column = 0; column < 4; column ++ ) results[idx].motif.count.push_back(baseCounts[idx][base]);
		}
	}
	vector<OutputMotif> motifs;
	selectRefinedOutputMotifs(&config, pipelines, results, motifs);
	require(motifs.size() == 2 && motifs[0].pipelineIdx == 1 && motifs[1].pipelineIdx == 2,
		"equal smoothed PWMs with different site counts must deduplicate in original Gibbs ranking");
	config.outputMotifNum = 1;
	selectRefinedOutputMotifs(&config, pipelines, results, motifs);
	require(motifs.size() == 1 && motifs[0].pipelineIdx == 1, "refined output limit follows deduplication and Gibbs ranking");
}

int main( void ) {
	checkFisherAndTies();
	checkBackgroundAndStrands();
	checkSubsetRefinement();
	checkRefinedDeduplication();
	cout << "All enrichment, Markov, strand, and subset-refinement tests passed.\n";
	return 0;
}
