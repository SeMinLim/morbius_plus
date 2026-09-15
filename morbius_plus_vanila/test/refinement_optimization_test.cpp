#include "Refinement.h"

#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;


static void require( bool condition, const char *message ) {
	if ( condition ) return;
	printf( "Refinement optimization test failed: %s\n", message );
	fflush( stdout );
	exit(1);
}

static void rebuildReferencePWM( const Config *config, const OutputMotif &motif,
				vector<double> &pwm ) {
	pwm.resize(4 * config->motifLength);
	double denominator = (double)motif.siteNum + 4 * PSEUDOCOUNT;
	for ( size_t idx = 0; idx < pwm.size(); idx ++ ) {
		pwm[idx] = ((double)motif.count[idx] + PSEUDOCOUNT) / denominator;
	}
}

// Serial, uncached control flow from commit 9745cd216a1537b887cccddd3a2ab4bf3bb9ba5e,
// using the current forward-only scan and refinement proposal limit.
// Always rescan the candidate, including an identical PWM; do not use the caches.
static void referenceRefineMotif( const Config *config, const Dataset *primary, const Dataset *control,
		const RefinementBackground *background, const OutputMotif &startingMotif,
		RefinementResult *result ) {
	OutputMotif initial = startingMotif;
	*result = RefinementResult{};
	result->motif = initial;
	result->motif.refined = true;
	std::vector<double> pwm;
	rebuildReferencePWM( config, result->motif, pwm );
	std::vector<RefinementSite> primarySites;
	std::vector<RefinementSite> controlSites;
	scanRefinementSites( config, primary, background, pwm, primarySites );
	scanRefinementSites( config, control, background, pwm, controlSites );
	result->selection = selectRefinementThreshold( primarySites, controlSites );
	if ( !result->selection.valid ) {
		result->motif.siteNum = 0;
		result->motif.count.assign( 4 * config->motifLength, 0 );
		result->motif.sitePresent.assign( primary->sequences.size(), 0 );
		result->motif.consensus.clear();
		result->terminationReason = "no_support";
		return;
	}
	result->initialLogPvalue = result->selection.logPvalue;
	for ( size_t iteration = 0; iteration < REFINEMENTMAXITERATIONS; iteration ++ ) {
		result->iterations++;
		OutputMotif candidate{};
		candidate.pipelineIdx = initial.pipelineIdx;
		candidate.refined = true;
		candidate.count.assign( 4 * config->motifLength, 0 );
		candidate.sitePresent.resize( primary->sequences.size() );
		candidate.offsets.resize( primary->sequences.size() );
		candidate.strands.resize( primary->sequences.size() );
		for ( size_t seqIdx = 0; seqIdx < primarySites.size(); seqIdx ++ ) {
			const RefinementSite &site = primarySites[seqIdx];
			bool present = site.score >= result->selection.scoreThreshold;
			candidate.sitePresent[seqIdx] = present ? 1 : 0;
			candidate.offsets[seqIdx] = site.offset;
			candidate.strands[seqIdx] = site.strand;
			if ( !present ) continue;
			candidate.siteNum++;
			for ( size_t column = 0; column < config->motifLength; column ++ ) {
				int symbol = getSiteSymbol( config, primary, primary->sequences[seqIdx],
					site.offset, site.strand, column );
				candidate.count[symbol * config->motifLength + column]++;
			}
		}
		if ( candidate.siteNum != result->selection.primarySiteNum ) {
			throw std::logic_error( "Refinement selected-site counts disagree with the threshold." );
		}
		candidate.consensus = buildConsensus( config, primary, candidate.count );
		std::vector<double> candidatePWM;
		rebuildReferencePWM( config, candidate, candidatePWM );
		std::vector<RefinementSite> candidatePrimarySites;
		std::vector<RefinementSite> candidateControlSites;
		scanRefinementSites( config, primary, background, candidatePWM, candidatePrimarySites );
		scanRefinementSites( config, control, background, candidatePWM, candidateControlSites );
		RefinementSelection candidateSelection = selectRefinementThreshold( candidatePrimarySites,
			candidateControlSites );
		if ( !candidateSelection.valid ) {
			result->terminationReason = "candidate_no_support";
			return;
		}
		// The returned PWM itself must improve enrichment, not just the sites used to fit it.
		if ( candidateSelection.logPvalue >= result->selection.logPvalue ) {
			result->terminationReason = "no_improvement";
			return;
		}
		result->fittingSelection = result->selection;
		result->scoringPWM = pwm;
		result->scoringSiteNum = result->motif.siteNum;
		result->selection = candidateSelection;
		result->motif = std::move( candidate );
		result->acceptedIterations++;
		result->acceptedLogPvalues.push_back( candidateSelection.logPvalue );
		pwm = std::move( candidatePWM );
		primarySites = std::move( candidatePrimarySites );
		controlSites = std::move( candidateControlSites );
	}
	result->terminationReason = "max_iterations";
}


static void compareSelection( const RefinementSelection &actual,
			     const RefinementSelection &expected ) {
	require(actual.valid == expected.valid, "selection validity");
	require(actual.scoreThreshold == expected.scoreThreshold, "exact score threshold");
	require(actual.logPvalue == expected.logPvalue, "exact Fisher log p-value");
	require(actual.primarySiteNum == expected.primarySiteNum, "Primary hit count");
	require(actual.controlSiteNum == expected.controlSiteNum, "Control hit count");
}

static void compareResult( const RefinementResult &actual, const RefinementResult &expected ) {
	compareSelection(actual.selection, expected.selection);
	compareSelection(actual.fittingSelection, expected.fittingSelection);
	require(actual.motif.pipelineIdx == expected.motif.pipelineIdx, "pipeline order");
	require(actual.motif.count == expected.motif.count, "exact PWM counts");
	require(actual.motif.consensus == expected.motif.consensus, "consensus");
	require(actual.motif.refined == expected.motif.refined, "refined flag");
	require(actual.motif.siteNum == expected.motif.siteNum, "fitting site count");
	require(actual.motif.sitePresent == expected.motif.sitePresent, "site presence mask");
	require(actual.motif.offsets == expected.motif.offsets, "all fitting offsets");
	require(actual.motif.strands == expected.motif.strands, "all fitting strands");
	for ( uint8_t strand : actual.motif.strands ) {
		require(strand == STRAND_FORWARD, "cached refinement returned a reverse site");
	}
	require(actual.initialLogPvalue == expected.initialLogPvalue, "initial log p-value");
	require(actual.iterations == expected.iterations, "proposal count");
	require(actual.acceptedIterations == expected.acceptedIterations, "accepted proposal count");
	require(actual.terminationReason == expected.terminationReason, "termination reason");
	require(actual.scoringPWM == expected.scoringPWM, "exact preceding scoring PWM");
	require(actual.scoringSiteNum == expected.scoringSiteNum, "preceding site count");
	require(actual.acceptedLogPvalues == expected.acceptedLogPvalues, "exact objective history");
}

static Dataset makeDataset( const vector<string> &sequences ) {
	Dataset dataset{};
	dataset.alphabet = "ACGT";
	dataset.alphabetSize = 4;
	dataset.symbolMap.assign(256, -1);
	for ( size_t idx = 0; idx < dataset.alphabet.size(); idx ++ ) {
		dataset.symbolMap[(unsigned char)dataset.alphabet[idx]] = (int)idx;
	}
	dataset.sequences = sequences;
	dataset.sequenceLength = sequences.empty() ? 0 : sequences[0].size();
	for ( size_t idx = 0; idx < sequences.size(); idx ++ ) {
		dataset.names.push_back("sequence_" + to_string(idx));
	}
	return dataset;
}

static OutputMotif makeInitial( const Config *config, const Dataset *primary,
			       const PipelineResult &pipeline, size_t pipelineIdx ) {
	OutputMotif initial{};
	initial.pipelineIdx = (int)pipelineIdx;
	initial.refined = true;
	initial.siteNum = primary->sequences.size();
	initial.sitePresent.assign(initial.siteNum, 1);
	initial.offsets = pipeline.bestOffsets;
	initial.strands = pipeline.bestStrands;
	buildResultCount(config, primary, initial.offsets, initial.strands, initial.count);
	initial.consensus = buildConsensus(config, primary, initial.count);
	return initial;
}

static void compareAllThreads( Config config, const Dataset &primary, const Dataset &control,
			      const vector<PipelineResult> &pipelines ) {
	RefinementBackground background;
	buildRefinementBackground(&control, &background);
	vector<RefinementResult> expected(pipelines.size());
	for ( size_t idx = 0; idx < pipelines.size(); idx ++ ) {
		OutputMotif initial = makeInitial(&config, &primary, pipelines[idx], idx);
		referenceRefineMotif(&config, &primary, &control, &background, initial, &expected[idx]);
		RefinementResult single;
		refineMotif(&config, &primary, &control, &background, initial, &single);
		compareResult(single, expected[idx]);
	}
	const int threadCounts[] = {1, 2, 4, 16, 0, 64};
	for ( size_t countIdx = 0; countIdx < sizeof(threadCounts) / sizeof(threadCounts[0]); countIdx ++ ) {
		config.threadNum = threadCounts[countIdx];
		RefinementBackground shared;
		vector<RefinementResult> actual;
		refinePipelineMotifs(&config, &primary, &control, pipelines, &shared, actual);
		require(shared.logProbability == background.logProbability, "shared Markov table");
		require(shared.letterFrequency == background.letterFrequency, "shared letter frequencies");
		require(actual.size() == expected.size(), "candidate count");
		for ( size_t idx = 0; idx < actual.size(); idx ++ ) compareResult(actual[idx], expected[idx]);
	}
}

static uint32_t fixtureWord( uint64_t *state ) {
	*state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
	return (uint32_t)(*state >> 32);
}

static string fixtureSequence( uint64_t *state, size_t length ) {
	string sequence(length, 'A');
	const string alphabet = "ACGT";
	for ( size_t idx = 0; idx < length; idx ++ ) sequence[idx] = alphabet[fixtureWord(state) % 4];
	return sequence;
}

static string reverseComplement( const string &sequence ) {
	const string alphabet = "ACGT";
	string reversed(sequence.size(), 'A');
	for ( size_t idx = 0; idx < sequence.size(); idx ++ ) {
		reversed[idx] = alphabet[3 - alphabet.find(sequence[sequence.size() - 1 - idx])];
	}
	return reversed;
}

static void checkGeneratedFixtures( void ) {
	const size_t widths[] = {1, 2, 4, 8, 10, 12, 16};
	for ( size_t caseIdx = 0; caseIdx < sizeof(widths) / sizeof(widths[0]); caseIdx ++ ) {
		Config config{};
		config.alphabetMode = ALPHABET_DNA;
		config.motifLength = widths[caseIdx];
		size_t length = caseIdx % 3 == 0 ? config.motifLength : 40;
		size_t primaryNum = 20 + caseIdx;
		size_t controlNum = 15 + 2 * caseIdx;
		uint64_t state = 1729 + caseIdx;
		vector<string> primarySequences;
		vector<string> controlSequences;
		string motif = fixtureSequence(&state, config.motifLength);
		vector<uint32_t> plantedOffsets(primaryNum);
		vector<uint8_t> plantedStrands(primaryNum);
		for ( size_t idx = 0; idx < primaryNum; idx ++ ) {
			// Exercise cache row boundaries with several internal sequence lengths.
			size_t sequenceLength = length + (caseIdx == 4 ? idx % 3 : 0);
			string sequence = fixtureSequence(&state, sequenceLength);
			plantedOffsets[idx] = fixtureWord(&state) % (sequenceLength - config.motifLength + 1);
			plantedStrands[idx] = (uint8_t)(idx % 2);
			if ( idx < 3 * primaryNum / 4 ) {
				sequence.replace(plantedOffsets[idx], config.motifLength,
					plantedStrands[idx] ? reverseComplement(motif) : motif);
			}
			primarySequences.push_back(sequence);
		}
		for ( size_t idx = 0; idx < controlNum; idx ++ ) {
			controlSequences.push_back(fixtureSequence(&state, length));
		}
		Dataset primary = makeDataset(primarySequences);
		Dataset control = makeDataset(controlSequences);
		vector<PipelineResult> pipelines(NUMPIPELINE);
		for ( size_t pipelineIdx = 0; pipelineIdx < pipelines.size(); pipelineIdx ++ ) {
			PipelineResult &pipeline = pipelines[pipelineIdx];
			pipeline.bestOffsets = plantedOffsets;
			pipeline.bestStrands.assign(primaryNum, STRAND_FORWARD);
			for ( size_t idx = 0; idx < primaryNum; idx ++ ) {
				if ( (idx + pipelineIdx) % 7 < pipelineIdx % 5 ) {
					pipeline.bestOffsets[idx] = fixtureWord(&state) %
						(primary.sequences[idx].size() - config.motifLength + 1);
					pipeline.bestStrands[idx] = STRAND_FORWARD;
				}
			}
		}
		compareAllThreads(config, primary, control, pipelines);
		if ( caseIdx == 3 ) compareAllThreads(config, primary, primary, pipelines);
	}
}

static void checkIdenticalPWMAndEmptyWork( void ) {
	Config config{};
	config.alphabetMode = ALPHABET_DNA;
	config.motifLength = 4;
	Dataset primary = makeDataset(vector<string>(8, "ACGAACGA"));
	Dataset control = makeDataset(vector<string>(5, "CCCCCCCC"));
	vector<PipelineResult> pipelines(NUMPIPELINE);
	for ( size_t idx = 0; idx < pipelines.size(); idx ++ ) {
		pipelines[idx].bestOffsets.assign(8, 0);
		pipelines[idx].bestStrands.assign(8, STRAND_FORWARD);
	}
	RefinementBackground background;
	buildRefinementBackground(&control, &background);
	OutputMotif initial = makeInitial(&config, &primary, pipelines[0], 0);
	RefinementResult reference;
	referenceRefineMotif(&config, &primary, &control, &background, initial, &reference);
	require(reference.selection.valid && reference.selection.primarySiteNum == primary.sequences.size(),
		"identical-PWM fixture has support at every original site");
	require(reference.iterations == 1 && reference.acceptedIterations == 0 &&
		reference.terminationReason == "no_improvement", "identical PWM must keep original diagnostics");
	compareAllThreads(config, primary, control, pipelines);
	pipelines.resize(3);
	compareAllThreads(config, primary, control, pipelines);
	pipelines.resize(1);
	compareAllThreads(config, primary, control, pipelines);
	pipelines.clear();
	compareAllThreads(config, primary, control, pipelines);
}

int main( void ) {
	checkGeneratedFixtures();
	checkIdenticalPWMAndEmptyWork();
	printf( "Refinement optimization tests passed: uncached equivalence, exact-PWM stop, and thread limits.\n" );
	return 0;
}
