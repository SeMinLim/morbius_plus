#include "MorbiusPlus.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <vector>
using namespace std;


#define PWLFRACTIONBITS 12
#define PWLCOEFFICIENTFRACTIONBITS 18
#define PWLRESIDUALBITS 8
#define PWLSEGMENTNUM 16


static const uint32_t LOG2PWLBASE[PWLSEGMENTNUM] = {
	0, 22928, 44545, 64993,
	84392, 102844, 120437, 137249,
	153344, 168783, 183616, 197889,
	211643, 224915, 237736, 250137
};

static const uint32_t LOG2PWLDELTA[PWLSEGMENTNUM] = {
	22928, 21617, 20448, 19399,
	18452, 17594, 16811, 16096,
	15439, 14833, 14273, 13754,
	13271, 12821, 12401, 12007
};

static const uint32_t EXP2PWLBASE[PWLSEGMENTNUM] = {
	262144, 251030, 240387, 230195,
	220436, 211090, 202141, 193571,
	185364, 177505, 169979, 162773,
	155872, 149263, 142935, 136875
};

static const uint32_t EXP2PWLDELTA[PWLSEGMENTNUM] = {
	11114, 10643, 10192, 9760,
	9346, 8950, 8570, 8207,
	7859, 7526, 7207, 6901,
	6608, 6328, 6060, 5803
};


// Calculate one 16-segment PWL value with the shared arithmetic structure
static uint32_t calculatePWLValue( const uint32_t *baseTable,
				   const uint32_t *deltaTable,
				   uint32_t fractionCode,
				   bool subtractCorrection ) {
	uint32_t segmentIdx = fractionCode >> PWLRESIDUALBITS;
	uint32_t residualMask = (1U << PWLRESIDUALBITS) - 1U;
	uint32_t residual = fractionCode & residualMask;
	uint64_t product = (uint64_t)deltaTable[segmentIdx] * (uint64_t)residual;
	uint32_t correction = (uint32_t)((product + (1U << (PWLRESIDUALBITS - 1))) >>
					 PWLRESIDUALBITS);

	if ( subtractCorrection ) return baseTable[segmentIdx] - correction;
	return baseTable[segmentIdx] + correction;
}

// Approximate log2(count) with the 16-segment LPM PWL datapath
static double calculateLog2PWL( uint32_t count ) {
	if ( count == 0 ) {
		printf( "Invalid zero BPM count for LPM conversion.\n" );
		fflush( stdout );
		exit(1);
	}

	uint32_t integerPart = 0;
	uint32_t normalizedCount = count;
	while ( normalizedCount > 1 ) {
		normalizedCount = normalizedCount >> 1;
		integerPart ++;
	}

	uint32_t leadingOne = 1U << integerPart;
	uint64_t fractionNumerator =
		(uint64_t)(count - leadingOne) << PWLFRACTIONBITS;
	uint32_t fractionCode = (uint32_t)(fractionNumerator / leadingOne);
	uint32_t fractionMaximum = (1U << PWLFRACTIONBITS) - 1U;
	if ( fractionCode > fractionMaximum ) fractionCode = fractionMaximum;

	uint32_t fractionQ18 = calculatePWLValue(LOG2PWLBASE,
						   LOG2PWLDELTA,
						   fractionCode,
						   false);
	uint32_t outputShift = PWLCOEFFICIENTFRACTIONBITS - PWLFRACTIONBITS;
	uint32_t fractionQ12 =
		(fractionQ18 + (1U << (outputShift - 1))) >> outputShift;
	if ( fractionQ12 >= (1U << PWLFRACTIONBITS) ) {
		integerPart ++;
		fractionQ12 = fractionQ12 - (1U << PWLFRACTIONBITS);
	}

	return (double)integerPart +
	       (double)fractionQ12 / (double)(1U << PWLFRACTIONBITS);
}

// Approximate 2^(-difference) with the 16-segment WeightProducer PWL datapath
static double calculateExp2NegativePWL( double difference ) {
	if ( difference < 0.0 ) {
		printf( "Invalid negative exponent difference for WeightProducer.\n" );
		fflush( stdout );
		exit(1);
	}

	uint64_t differenceQ12 =
		(uint64_t)llround(difference * (double)(1U << PWLFRACTIONBITS));
	uint32_t integerPart = (uint32_t)(differenceQ12 >> PWLFRACTIONBITS);
	uint32_t fractionCode =
		(uint32_t)(differenceQ12 & ((1U << PWLFRACTIONBITS) - 1U));
	uint32_t fractionQ18 = calculatePWLValue(EXP2PWLBASE,
						   EXP2PWLDELTA,
						   fractionCode,
						   true);

	if ( integerPart > PWLCOEFFICIENTFRACTIONBITS ) return 0.0;
	uint32_t weightQ18 = fractionQ18 >> integerPart;
	return (double)weightQ18 /
	       (double)(1U << PWLCOEFFICIENTFRACTIONBITS);
}

// Initialize offsets with fixed-rate support-guided initialization
void initializeOffsets( const Config *config,
			const Dataset *dataset,
			const SeedModel *seedModel,
			int pipelineIdx,
			vector<uint32_t> &offsets ) {
	RandomGenerator initializationRandom;
	initializeRandomGenerator(&initializationRandom,
				  config->randomSeed +
				  (uint64_t)(pipelineIdx + 1) * 0x9e3779b97f4a7c15ULL);

	offsets.resize(dataset->sequences.size());
	uint32_t candidateNum = (uint32_t)(dataset->sequenceLength - config->motifLength + 1);
	for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx ++ ) {
		uint32_t randomOffset = randomBounded(&initializationRandom, candidateNum);
		uint32_t guidedOffset = randomOffset;
		bool guidedValid = false;

		if ( seedModel->valid && seedModel->guidedOffsets[seqIdx].empty() == false ) {
			uint32_t guidedIdx = randomBounded(&initializationRandom,
							   (uint32_t)seedModel->guidedOffsets[seqIdx].size());
			guidedOffset = seedModel->guidedOffsets[seqIdx][guidedIdx];
			guidedValid = true;
		}

		uint32_t guidanceBit = randomWord(&initializationRandom) & 0x1U;
		if ( guidanceBit == 1 && guidedValid ) offsets[seqIdx] = guidedOffset;
		else offsets[seqIdx] = randomOffset;
	}
}

// Build a complete BPM
void buildBPM( const Config *config,
	       const Dataset *dataset,
	       const vector<uint32_t> &offsets,
	       vector<uint32_t> &bpm ) {
	bpm.assign((size_t)dataset->alphabetSize * config->motifLength, PSEUDOCOUNT);
	for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx ++ ) {
		for ( size_t column = 0; column < config->motifLength; column ++ ) {
			int symbol = dataset->symbolMap[(unsigned char)dataset->sequences[seqIdx][offsets[seqIdx] + column]];
			bpm[(size_t)symbol * config->motifLength + column] ++;
		}
	}
}

// Build a complete LPM
void buildLPM( const vector<uint32_t> &bpm, vector<double> &lpm ) {
	lpm.resize(bpm.size());
	for ( size_t i = 0; i < bpm.size(); i ++ ) {
		lpm[i] = calculateLog2PWL(bpm[i]);
	}
}

// Calculate the overall consensus agreement score
uint64_t calculateAgreementScore( const Config *config,
				  const Dataset *dataset,
				  const vector<uint32_t> &bpm ) {
	uint64_t score = 0;
	for ( size_t column = 0; column < config->motifLength; column ++ ) {
		uint32_t columnMaximum = 0;
		for ( int symbol = 0; symbol < dataset->alphabetSize; symbol ++ ) {
			columnMaximum = max(columnMaximum,
					    bpm[(size_t)symbol * config->motifLength + column]);
		}
		score = score + columnMaximum;
	}
	return score;
}

// Calculate a normalized agreement score
double calculateNormalizedScore( const Config *config,
					 const Dataset *dataset,
					 uint64_t score ) {
	double minimumPseudocountScore = (double)config->motifLength;
	double denominator = (double)dataset->sequences.size() * (double)config->motifLength;
	return ((double)score - minimumPseudocountScore) / denominator;
}

// Calculate the raw score threshold
uint64_t calculateRawScoreThreshold( const Config *config, const Dataset *dataset ) {
	double score = (double)config->motifLength +
		       config->scoreThreshold *
		       (double)dataset->sequences.size() *
		       (double)config->motifLength;
	return (uint64_t)ceil(score);
}

// Calculate one candidate LogProb
double calculateCandidateLogProb( const Config *config,
					  const Dataset *dataset,
					  const string &sequence,
					  size_t candidateOffset,
					  const vector<double> &lpm ) {
	double logProb = 0.0;
	for ( size_t column = 0; column < config->motifLength; column ++ ) {
		int symbol = dataset->symbolMap[(unsigned char)sequence[candidateOffset + column]];
		logProb = logProb + lpm[(size_t)symbol * config->motifLength + column];
	}
	return logProb;
}

// Sample a candidate with hierarchical inverse-CDF
uint32_t sampleCandidate( const Config *config,
			  const Dataset *dataset,
			  const string &sequence,
			  const vector<double> &lpm,
			  RandomGenerator *randomGenerator ) {
	size_t candidateNum = sequence.size() - config->motifLength + 1;
	vector<SegmentSummary> segmentSummaries;
	segmentSummaries.reserve((candidateNum + SEGMENTSIZE - 1) / SEGMENTSIZE);
	vector<double> logProb(SEGMENTSIZE, 0.0);
	vector<double> weight(SEGMENTSIZE, 0.0);

	bool globalInitialized = false;
	int globalExponent = 0;
	double globalMass = 0.0;

	for ( size_t segmentStart = 0; segmentStart < candidateNum; segmentStart += SEGMENTSIZE ) {
		size_t segmentCandidateNum = min((size_t)SEGMENTSIZE, candidateNum - segmentStart);
		double segmentMaximum = -numeric_limits<double>::infinity();

		for ( size_t localIdx = 0; localIdx < segmentCandidateNum; localIdx ++ ) {
			logProb[localIdx] = calculateCandidateLogProb(config,
									  dataset,
									  sequence,
									  segmentStart + localIdx,
									  lpm);
			segmentMaximum = max(segmentMaximum, logProb[localIdx]);
		}

		int segmentExponent = (int)ceil(segmentMaximum);
		double segmentMass = 0.0;
		for ( size_t localIdx = 0; localIdx < segmentCandidateNum; localIdx ++ ) {
			double difference = (double)segmentExponent - logProb[localIdx];
			weight[localIdx] = calculateExp2NegativePWL(difference);
			segmentMass = segmentMass + weight[localIdx];
		}

		double localThreshold = randomUnit(randomGenerator) * segmentMass;
		double localCumulative = 0.0;
		size_t selectedLocalOffset = segmentCandidateNum - 1;
		for ( size_t localIdx = 0; localIdx < segmentCandidateNum; localIdx ++ ) {
			localCumulative = localCumulative + weight[localIdx];
			if ( localCumulative > localThreshold ) {
				selectedLocalOffset = localIdx;
				break;
			}
		}

		SegmentSummary summary;
		summary.startOffset = segmentStart;
		summary.exponent = segmentExponent;
		summary.mass = segmentMass;
		summary.localOffset = selectedLocalOffset;
		segmentSummaries.push_back(summary);

		if ( globalInitialized == false ) {
			globalExponent = segmentExponent;
			globalMass = segmentMass;
			globalInitialized = true;
		} else if ( segmentExponent <= globalExponent ) {
			globalMass = globalMass + ldexp(segmentMass, segmentExponent - globalExponent);
		} else {
			globalMass = ldexp(globalMass, globalExponent - segmentExponent) + segmentMass;
			globalExponent = segmentExponent;
		}
	}

	double globalThreshold = randomUnit(randomGenerator) * globalMass;
	double globalCumulative = 0.0;
	for ( size_t segmentIdx = 0; segmentIdx < segmentSummaries.size(); segmentIdx ++ ) {
		const SegmentSummary &summary = segmentSummaries[segmentIdx];
		double alignedMass = ldexp(summary.mass, summary.exponent - globalExponent);
		globalCumulative = globalCumulative + alignedMass;
		if ( globalCumulative > globalThreshold || segmentIdx + 1 == segmentSummaries.size() ) {
			return (uint32_t)(summary.startOffset + summary.localOffset);
		}
	}

	return 0;
}

// Prepare the initial leave-one-out BPM and LPM
void prepareInitialPipelineState( const Config *config,
				  const Dataset *dataset,
				  PipelineState *state ) {
	buildBPM(config, dataset, state->offsets, state->bpm);
	state->currentScore = calculateAgreementScore(config, dataset, state->bpm);
	state->bestScore = state->currentScore;
	state->bestOffsets = state->offsets;
	state->updateNum = 0;
	state->thresholdReached = false;

	// Remove sequence 0 to form the first leave-one-out BPM
	for ( size_t column = 0; column < config->motifLength; column ++ ) {
		int symbol = dataset->symbolMap[(unsigned char)dataset->sequences[0][state->offsets[0] + column]];
		state->bpm[(size_t)symbol * config->motifLength + column] --;
	}
	buildLPM(state->bpm, state->lpm);
}

// Run one Gibbs pipeline
void runPipeline( const Config *config,
		  const Dataset *dataset,
		  const SeedModel *seedModel,
		  int pipelineIdx,
		  uint64_t rawScoreThreshold,
		  PipelineResult *result ) {
	double startTime = timeChecker();
	PipelineState state;
	initializeOffsets(config, dataset, seedModel, pipelineIdx, state.offsets);
	initializeRandomGenerator(&state.randomGenerator,
				  config->randomSeed ^
				  ((uint64_t)(pipelineIdx + 1) * 0xd2b74407b1ce6e93ULL));
	prepareInitialPipelineState(config, dataset, &state);

	if ( state.bestScore >= rawScoreThreshold ) state.thresholdReached = true;
	uint64_t maxUpdateNum = config->maxUpdateNum;
	for ( size_t seqIdx = 0;
	      state.thresholdReached == false && state.updateNum < maxUpdateNum;
	      seqIdx = (seqIdx + 1) % dataset->sequences.size() ) {
		uint32_t newOffset = sampleCandidate(config,
						     dataset,
						     dataset->sequences[seqIdx],
						     state.lpm,
						     &state.randomGenerator);

		// Add the new tentative motif to create the complete BPM
		for ( size_t column = 0; column < config->motifLength; column ++ ) {
			int symbol = dataset->symbolMap[(unsigned char)dataset->sequences[seqIdx][newOffset + column]];
			state.bpm[(size_t)symbol * config->motifLength + column] ++;
		}
		state.offsets[seqIdx] = newOffset;
		state.currentScore = calculateAgreementScore(config, dataset, state.bpm);
		state.updateNum ++;

		// Pipeline-local ScoreComparator
		if ( state.currentScore > state.bestScore ) {
			state.bestScore = state.currentScore;
			state.bestOffsets = state.offsets;
		}
		if ( state.bestScore >= rawScoreThreshold ) {
			state.thresholdReached = true;
			break;
		}
		if ( state.updateNum >= maxUpdateNum ) break;

		// Remove the next old tentative motif and selectively update the LPM
		size_t nextSeqIdx = (seqIdx + 1) % dataset->sequences.size();
		uint32_t nextOldOffset = state.offsets[nextSeqIdx];
		for ( size_t column = 0; column < config->motifLength; column ++ ) {
			int newSymbol = dataset->symbolMap[(unsigned char)dataset->sequences[seqIdx][newOffset + column]];
			int nextOldSymbol = dataset->symbolMap[(unsigned char)dataset->sequences[nextSeqIdx][nextOldOffset + column]];

			state.bpm[(size_t)nextOldSymbol * config->motifLength + column] --;
			if ( newSymbol != nextOldSymbol ) {
				state.lpm[(size_t)newSymbol * config->motifLength + column] =
					calculateLog2PWL(state.bpm[(size_t)newSymbol * config->motifLength + column]);
				state.lpm[(size_t)nextOldSymbol * config->motifLength + column] =
					calculateLog2PWL(state.bpm[(size_t)nextOldSymbol * config->motifLength + column]);
			}
		}
	}

	result->bestOffsets = state.bestOffsets;
	result->bestScore = state.bestScore;
	result->updateNum = state.updateNum;
	result->thresholdReached = state.thresholdReached;
	result->elapsedTime = timeChecker() - startTime;
}

// Run all Gibbs pipelines
void runPipelines( const Config *config,
		   const Dataset *dataset,
		   const SeedModel *seedModel,
		   vector<PipelineResult> &pipelineResults ) {
	pipelineResults.resize((size_t)NUMPIPELINE);
	uint64_t rawScoreThreshold = calculateRawScoreThreshold(config, dataset);
	int threadNum = config->threadNum;
	if ( threadNum == 0 ) threadNum = NUMPIPELINE;
	threadNum = min(threadNum, NUMPIPELINE);
	if ( threadNum < 1 ) threadNum = 1;

	for ( int batchStart = 0; batchStart < NUMPIPELINE; batchStart += threadNum ) {
		int batchNum = min(threadNum, NUMPIPELINE - batchStart);
		vector<thread> workers;
		workers.reserve((size_t)batchNum);
		for ( int workerIdx = 0; workerIdx < batchNum; workerIdx ++ ) {
			int pipelineIdx = batchStart + workerIdx;
			workers.emplace_back(runPipeline,
					     config,
					     dataset,
					     seedModel,
					     pipelineIdx,
					     rawScoreThreshold,
					     &pipelineResults[pipelineIdx]);
		}
		for ( size_t workerIdx = 0; workerIdx < workers.size(); workerIdx ++ ) workers[workerIdx].join();
		for ( int workerIdx = 0; workerIdx < batchNum; workerIdx ++ ) {
			int pipelineIdx = batchStart + workerIdx;
			printf( "[STEP 3] Pipeline %d is done. Score: %lu, Updates: %lu\n",
				pipelineIdx,
				(unsigned long)pipelineResults[pipelineIdx].bestScore,
				(unsigned long)pipelineResults[pipelineIdx].updateNum
			);
		}
		fflush( stdout );
	}
}

// Select the pipeline with the highest best score
int selectBestPipeline( const vector<PipelineResult> &pipelineResults ) {
	int bestPipelineIdx = 0;
	for ( int pipelineIdx = 1; pipelineIdx < (int)pipelineResults.size(); pipelineIdx ++ ) {
		if ( pipelineResults[pipelineIdx].bestScore > pipelineResults[bestPipelineIdx].bestScore ) {
			bestPipelineIdx = pipelineIdx;
		}
	}
	return bestPipelineIdx;
}
