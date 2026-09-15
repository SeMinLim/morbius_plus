#include "Refinement.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <thread>


namespace {

struct FisherTable {
	std::vector<long double> logFactorial;

	explicit FisherTable( size_t total ) : logFactorial( total + 1, 0.0L ) {
		// Compensated accumulation keeps cancellation in log binomial terms small.
		long double sum = 0.0L;
		long double correction = 0.0L;
		for ( size_t idx = 2; idx <= total; idx++ ) {
			long double term = std::log( (long double)idx ) - correction;
			long double next = sum + term;
			correction = (next - sum) - term;
			sum = next;
			logFactorial[idx] = sum;
		}
	}

	long double logChoose( size_t total, size_t chosen ) const {
		return logFactorial[total] - logFactorial[chosen] - logFactorial[total - chosen];
	}

	long double logMass( size_t hits, size_t totalHits, size_t primaryNum, size_t total ) const {
		return logChoose( totalHits, hits ) + logChoose( total - totalHits, primaryNum - hits )
			- logChoose( total, primaryNum );
	}
};

long double upperTailLogPvalue( const FisherTable &table, size_t hits,
		size_t totalHits, size_t primaryNum, size_t total ) {
	size_t last = std::min( totalHits, primaryNum );
	long double term = 1.0L;
	long double sum = 1.0L;
	const long double tolerance = 8.0L * std::numeric_limits<long double>::epsilon();
	for ( size_t idx = hits; idx < last; idx++ ) {
		long double ratio = ((long double)(totalHits - idx) * (primaryNum - idx)) /
			((long double)(idx + 1) * (total - totalHits - (primaryNum - idx) + 1));
		// Above the mean these ratios decrease. This bounds the entire omitted tail.
		if ( ratio < 1.0L && term * ratio / (1.0L - ratio) <= tolerance * sum ) break;
		term *= ratio;
		sum += term;
	}
	return std::min( 0.0L, table.logMass( hits, totalHits, primaryNum, total ) + std::log( sum ) );
}

long double lowerTailLogPvalue( const FisherTable &table, size_t hits,
		size_t totalHits, size_t primaryNum, size_t total ) {
	size_t first = totalHits > total - primaryNum ? totalHits - (total - primaryNum) : 0;
	long double term = 1.0L;
	long double sum = 1.0L;
	const long double tolerance = 8.0L * std::numeric_limits<long double>::epsilon();
	for ( size_t idx = hits; idx > first; idx-- ) {
		long double ratio = ((long double)idx * (total - totalHits - (primaryNum - idx))) /
			((long double)(totalHits - idx + 1) * (primaryNum - idx + 1));
		if ( ratio < 1.0L && term * ratio / (1.0L - ratio) <= tolerance * sum ) break;
		term *= ratio;
		sum += term;
	}
	return std::min( 0.0L, table.logMass( hits, totalHits, primaryNum, total ) + std::log( sum ) );
}

long double evaluateFisherLogPvalue( const FisherTable &table, size_t primaryHits,
		size_t controlHits, size_t primaryNum, size_t controlNum ) {
	size_t total = primaryNum + controlNum;
	size_t totalHits = primaryHits + controlHits;
	size_t first = totalHits > controlNum ? totalHits - controlNum : 0;
	if ( primaryHits <= first || total == 0 ) return 0.0L;
	if ( (long double)primaryHits * total > (long double)totalHits * primaryNum ) {
		return upperTailLogPvalue( table, primaryHits, totalHits, primaryNum, total );
	}
	// Below the mean, sum the shorter lower tail and take its complement stably.
	long double lowerLogPvalue = lowerTailLogPvalue( table, primaryHits - 1,
		totalHits, primaryNum, total );
	return std::log( -std::expm1( lowerLogPvalue ) );
}

void validatePopulationSizes( size_t primaryNum, size_t controlNum ) {
	if ( controlNum == std::numeric_limits<size_t>::max() ||
		primaryNum > std::numeric_limits<size_t>::max() - controlNum - 1 ) {
		throw std::overflow_error( "Refinement population size is too large." );
	}
}

// One shared, immutable background score per sequence, strand, and offset.
struct RefinementBackgroundCache {
	std::vector<size_t> sequenceStart;
	std::vector<double> logProbability;
};

struct ScoredSequence {
	double score;
	bool primary;
};

void rebuildRefinementPWM( const Config *config, const OutputMotif &motif,
		std::vector<double> &pwm ) {
	pwm.resize( 4 * config->motifLength );
	double denominator = (double)motif.siteNum + 4 * PSEUDOCOUNT;
	for ( size_t idx = 0; idx < pwm.size(); idx++ ) {
		pwm[idx] = ((double)motif.count[idx] + PSEUDOCOUNT) / denominator;
	}
}

}


void buildRefinementBackground( const Dataset *control, RefinementBackground *background ) {
	if ( control->alphabetSize != 4 || control->sequences.empty() ) {
		throw std::invalid_argument( "Refinement requires a nonempty DNA Control dataset." );
	}
	background->logProbability.resize( REFINEMENTBACKGROUNDORDER + 1 );
	for ( size_t order = 0; order <= REFINEMENTBACKGROUNDORDER; order++ ) {
		size_t tupleLength = order + 1;
		size_t tupleNum = (size_t)1 << (2 * tupleLength);
		std::vector<long double> count( tupleNum, 1.0L / tupleNum );
		for ( size_t seqIdx = 0; seqIdx < control->sequences.size(); seqIdx++ ) {
			const std::string &sequence = control->sequences[seqIdx];
			if ( sequence.size() < tupleLength ) continue;
			for ( size_t offset = 0; offset <= sequence.size() - tupleLength; offset++ ) {
				uint32_t forward = 0;
				uint32_t reverse = 0;
				for ( size_t column = 0; column < tupleLength; column++ ) {
					forward = forward * 4 + control->symbolMap[(unsigned char)sequence[offset + column]];
					reverse = reverse * 4 + 3 - control->symbolMap[(unsigned char)sequence[offset + tupleLength - 1 - column]];
				}
				count[forward] += 1.0L;
				count[reverse] += 1.0L;
			}
		}
		background->logProbability[order].resize( tupleNum );
		for ( size_t context = 0; context < tupleNum / 4; context++ ) {
			long double total = 0.0L;
			for ( size_t base = 0; base < 4; base++ ) total += count[context * 4 + base];
			for ( size_t base = 0; base < 4; base++ ) {
				background->logProbability[order][context * 4 + base] =
					(double)std::log2( count[context * 4 + base] / total );
			}
		}
	}
	background->letterFrequency.resize( 4 );
	for ( size_t base = 0; base < 4; base++ ) {
		background->letterFrequency[base] = std::exp2( background->logProbability[0][base] );
	}
}


double calculateRefinementBackgroundLogProbability( const Config *config,
		const Dataset *dataset, const RefinementBackground *background,
		const std::string &sequence, uint32_t offset, uint8_t strand ) {
	uint32_t context = 0;
	double logProbability = 0.0;
	for ( size_t column = 0; column < config->motifLength; column++ ) {
		int symbol = getSiteSymbol( config, dataset, sequence, offset, strand, column );
		size_t order = std::min( column, (size_t)REFINEMENTBACKGROUNDORDER );
		logProbability += background->logProbability[order][context * 4 + symbol];
		context = (context * 4 + symbol) & 15;
	}
	return logProbability;
}


void scanRefinementSites( const Config *config, const Dataset *dataset,
		const RefinementBackground *background, const std::vector<double> &pwm,
		std::vector<RefinementSite> &sites ) {
	if ( config->alphabetMode != ALPHABET_DNA || dataset->alphabetSize != 4 ||
		config->motifLength == 0 || config->motifLength > dataset->sequenceLength ||
		pwm.size() != 4 * config->motifLength ) {
		throw std::invalid_argument( "Invalid DNA PWM or sequence length for refinement." );
	}
	std::vector<double> logPWM( pwm.size() );
	for ( size_t idx = 0; idx < pwm.size(); idx++ ) {
		if ( !std::isfinite( pwm[idx] ) || pwm[idx] < 0.0 ) {
			throw std::invalid_argument( "Refinement PWM probabilities must be finite and nonnegative." );
		}
		logPWM[idx] = pwm[idx] == 0.0 ? -std::numeric_limits<double>::infinity() : std::log2( pwm[idx] );
	}
	sites.resize( dataset->sequences.size() );
	for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx++ ) {
		const std::string &sequence = dataset->sequences[seqIdx];
		if ( sequence.size() < config->motifLength || sequence.size() - config->motifLength > UINT32_MAX ) {
			throw std::invalid_argument( "Invalid sequence length for refinement site offsets." );
		}
		RefinementSite best = { -std::numeric_limits<double>::infinity(), 0, STRAND_FORWARD };
		// Iteration order resolves exact ties: forward strand, then the lowest offset.
		for ( uint8_t strand = STRAND_FORWARD; strand <= STRAND_REVERSE; strand++ ) {
			for ( size_t offset = 0; offset <= sequence.size() - config->motifLength; offset++ ) {
				uint32_t context = 0;
				double motifLogProbability = 0.0;
				double backgroundLogProbability = 0.0;
				for ( size_t column = 0; column < config->motifLength; column++ ) {
					int symbol = getSiteSymbol( config, dataset, sequence, offset, strand, column );
					size_t order = std::min( column, (size_t)REFINEMENTBACKGROUNDORDER );
					motifLogProbability += logPWM[symbol * config->motifLength + column];
					backgroundLogProbability += background->logProbability[order][context * 4 + symbol];
					context = (context * 4 + symbol) & 15;
				}
				double score = motifLogProbability - backgroundLogProbability;
				if ( score > best.score ) {
					best.score = score;
					best.offset = (uint32_t)offset;
					best.strand = strand;
				}
			}
		}
		sites[seqIdx] = best;
	}
}


namespace {

void buildRefinementBackgroundCache( const Config *config, const Dataset *dataset,
		const RefinementBackground *background, RefinementBackgroundCache *cache ) {
	if ( config->alphabetMode != ALPHABET_DNA || dataset->alphabetSize != 4 ||
		config->motifLength == 0 || config->motifLength > dataset->sequenceLength ) {
		throw std::invalid_argument( "Invalid DNA sequence length for refinement cache." );
	}
	if ( dataset->sequences.size() >= cache->sequenceStart.max_size() ) {
		throw std::length_error( "Too many sequences for the refinement cache." );
	}
	cache->sequenceStart.resize( dataset->sequences.size() + 1 );
	size_t scoreNum = 0;
	for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx ++ ) {
		size_t length = dataset->sequences[seqIdx].size();
		if ( length < config->motifLength || length - config->motifLength > UINT32_MAX ) {
			throw std::invalid_argument( "Invalid sequence length for refinement site offsets." );
		}
		size_t offsetNum = length - config->motifLength + 1;
		if ( offsetNum > (cache->logProbability.max_size() - scoreNum) / 2 ) {
			throw std::length_error( "Refinement background cache is too large." );
		}
		cache->sequenceStart[seqIdx] = scoreNum;
		scoreNum += 2 * offsetNum;
	}
	cache->sequenceStart[dataset->sequences.size()] = scoreNum;
	cache->logProbability.resize( scoreNum );
	for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx ++ ) {
		const std::string &sequence = dataset->sequences[seqIdx];
		size_t scoreIdx = cache->sequenceStart[seqIdx];
		for ( uint8_t strand = STRAND_FORWARD; strand <= STRAND_REVERSE; strand ++ ) {
			for ( size_t offset = 0; offset <= sequence.size() - config->motifLength; offset ++ ) {
				cache->logProbability[scoreIdx] = calculateRefinementBackgroundLogProbability(
					config, dataset, background, sequence, (uint32_t)offset, strand );
				scoreIdx ++;
			}
		}
	}
}

void scanCachedRefinementSites( const Config *config, const Dataset *dataset,
		const RefinementBackgroundCache *cache, const std::vector<double> &pwm,
		std::vector<RefinementSite> &sites ) {
	if ( pwm.size() != 4 * config->motifLength ) {
		throw std::invalid_argument( "Invalid DNA PWM or sequence length for refinement." );
	}
	std::vector<double> logPWM( pwm.size() );
	for ( size_t idx = 0; idx < pwm.size(); idx ++ ) {
		if ( !std::isfinite( pwm[idx] ) || pwm[idx] < 0.0 ) {
			throw std::invalid_argument( "Refinement PWM probabilities must be finite and nonnegative." );
		}
		logPWM[idx] = pwm[idx] == 0.0 ? -std::numeric_limits<double>::infinity() : std::log2( pwm[idx] );
	}
	sites.resize( dataset->sequences.size() );
	for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx ++ ) {
		const std::string &sequence = dataset->sequences[seqIdx];
		size_t scoreIdx = cache->sequenceStart[seqIdx];
		RefinementSite best = { -std::numeric_limits<double>::infinity(), 0, STRAND_FORWARD };
		// Keep the original sum order and ties: forward strand, then lowest offset.
		for ( uint8_t strand = STRAND_FORWARD; strand <= STRAND_REVERSE; strand ++ ) {
			for ( size_t offset = 0; offset <= sequence.size() - config->motifLength; offset ++ ) {
				double motifLogProbability = 0.0;
				for ( size_t column = 0; column < config->motifLength; column ++ ) {
					int symbol = getSiteSymbol( config, dataset, sequence, offset, strand, column );
					motifLogProbability += logPWM[symbol * config->motifLength + column];
				}
				double score = motifLogProbability - cache->logProbability[scoreIdx];
				scoreIdx ++;
				if ( score > best.score ) {
					best.score = score;
					best.offset = (uint32_t)offset;
					best.strand = strand;
				}
			}
		}
		sites[seqIdx] = best;
	}
}

}


double fisherLogPvalue( size_t primaryHits, size_t controlHits,
		size_t primaryNum, size_t controlNum ) {
	if ( primaryHits > primaryNum || controlHits > controlNum ) {
		throw std::invalid_argument( "Fisher hit counts exceed their sequence populations." );
	}
	validatePopulationSizes( primaryNum, controlNum );
	FisherTable table( primaryNum + controlNum );
	return (double)evaluateFisherLogPvalue( table, primaryHits, controlHits, primaryNum, controlNum );
}


namespace {

RefinementSelection selectRefinementThresholdWithTable( const std::vector<RefinementSite> &primarySites,
		const std::vector<RefinementSite> &controlSites, const FisherTable &table ) {
	RefinementSelection best;
	if ( primarySites.empty() || controlSites.empty() ) return best;
	validatePopulationSizes( primarySites.size(), controlSites.size() );
	size_t total = primarySites.size() + controlSites.size();
	std::vector<ScoredSequence> scores;
	scores.reserve( total );
	for ( size_t idx = 0; idx < primarySites.size(); idx++ ) {
		if ( std::isfinite( primarySites[idx].score ) && primarySites[idx].score > 0.0 ) {
			scores.push_back( { primarySites[idx].score, true } );
		}
	}
	for ( size_t idx = 0; idx < controlSites.size(); idx++ ) {
		if ( std::isfinite( controlSites[idx].score ) && controlSites[idx].score > 0.0 ) {
			scores.push_back( { controlSites[idx].score, false } );
		}
	}
	std::sort( scores.begin(), scores.end(), []( const ScoredSequence &left, const ScoredSequence &right ) {
		return left.score > right.score;
	} );
	size_t primaryHits = 0;
	size_t controlHits = 0;
	long double bestLogPvalue = 0.0L;
	for ( size_t idx = 0; idx < scores.size(); ) {
		double threshold = scores[idx].score;
		do {
			if ( scores[idx].primary ) primaryHits++;
			else controlHits++;
			idx++;
		} while ( idx < scores.size() && scores[idx].score == threshold );
		// Never split a score tie or retain a cut depleted in Primary.
		if ( (long double)primaryHits * controlSites.size() <=
			(long double)controlHits * primarySites.size() ) continue;
		long double logMass = table.logMass( primaryHits, primaryHits + controlHits,
			primarySites.size(), total );
		// A tail is at least its first mass, so this cut cannot improve the current best.
		if ( logMass >= bestLogPvalue ) continue;
		long double logPvalue = evaluateFisherLogPvalue( table, primaryHits, controlHits,
			primarySites.size(), controlSites.size() );
		if ( logPvalue < bestLogPvalue ) {
			best.valid = true;
			best.scoreThreshold = threshold;
			best.logPvalue = (double)logPvalue;
			best.primarySiteNum = primaryHits;
			best.controlSiteNum = controlHits;
			bestLogPvalue = logPvalue;
		}
	}
	return best;
}


void refineMotifWithCache( const Config *config, const Dataset *primary, const Dataset *control,
		const RefinementBackgroundCache *primaryCache, const RefinementBackgroundCache *controlCache,
		const FisherTable &table, const OutputMotif &startingMotif, RefinementResult *result ) {
	OutputMotif initial = startingMotif;
	*result = RefinementResult{};
	result->motif = initial;
	result->motif.refined = true;
	std::vector<double> pwm;
	rebuildRefinementPWM( config, result->motif, pwm );
	std::vector<RefinementSite> primarySites;
	std::vector<RefinementSite> controlSites;
	scanCachedRefinementSites( config, primary, primaryCache, pwm, primarySites );
	scanCachedRefinementSites( config, control, controlCache, pwm, controlSites );
	result->selection = selectRefinementThresholdWithTable( primarySites, controlSites, table );
	if ( !result->selection.valid ) {
		result->motif.siteNum = 0;
		result->motif.count.assign( 4 * config->motifLength, 0 );
		result->motif.sitePresent.assign( primary->sequences.size(), 0 );
		result->motif.consensus.clear();
		result->terminationReason = "no_support";
		return;
	}
	result->initialLogPvalue = result->selection.logPvalue;
	for ( size_t iteration = 0; iteration < REFINEMENTMAXITERATIONS; iteration++ ) {
		result->iterations++;
		OutputMotif candidate{};
		candidate.pipelineIdx = initial.pipelineIdx;
		candidate.refined = true;
		candidate.count.assign( 4 * config->motifLength, 0 );
		candidate.sitePresent.resize( primary->sequences.size() );
		candidate.offsets.resize( primary->sequences.size() );
		candidate.strands.resize( primary->sequences.size() );
		for ( size_t seqIdx = 0; seqIdx < primarySites.size(); seqIdx++ ) {
			const RefinementSite &site = primarySites[seqIdx];
			bool present = site.score >= result->selection.scoreThreshold;
			candidate.sitePresent[seqIdx] = present ? 1 : 0;
			candidate.offsets[seqIdx] = site.offset;
			candidate.strands[seqIdx] = site.strand;
			if ( !present ) continue;
			candidate.siteNum++;
			for ( size_t column = 0; column < config->motifLength; column++ ) {
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
		rebuildRefinementPWM( config, candidate, candidatePWM );
		// Identical PWM gives the same scan and Fisher result; preserve the old exit status.
		if ( candidatePWM == pwm ) {
			result->terminationReason = "no_improvement";
			return;
		}
		std::vector<RefinementSite> candidatePrimarySites;
		std::vector<RefinementSite> candidateControlSites;
		scanCachedRefinementSites( config, primary, primaryCache, candidatePWM, candidatePrimarySites );
		scanCachedRefinementSites( config, control, controlCache, candidatePWM, candidateControlSites );
		RefinementSelection candidateSelection = selectRefinementThresholdWithTable( candidatePrimarySites,
			candidateControlSites, table );
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


struct RefinementWork {
	const Config *config;
	const Dataset *primary;
	const Dataset *control;
	const std::vector<PipelineResult> *pipelineResults;
	const RefinementBackgroundCache *primaryCache;
	const RefinementBackgroundCache *controlCache;
	const FisherTable *table;
	std::vector<RefinementResult> *results;
	std::vector<std::exception_ptr> *errors;
	std::atomic<size_t> nextPipelineIdx{0};
};

// Each worker owns one candidate at a time; all shared lookup tables are read-only.
void refinePipelineWorker( RefinementWork *work ) {
	for ( ;; ) {
		size_t pipelineIdx = work->nextPipelineIdx.fetch_add( 1, std::memory_order_relaxed );
		if ( pipelineIdx >= work->pipelineResults->size() ) return;
		try {
			const PipelineResult &pipeline = (*work->pipelineResults)[pipelineIdx];
			OutputMotif initial{};
			initial.pipelineIdx = (int)pipelineIdx;
			initial.refined = true;
			initial.siteNum = work->primary->sequences.size();
			initial.sitePresent.assign( work->primary->sequences.size(), 1 );
			initial.offsets = pipeline.bestOffsets;
			initial.strands = pipeline.bestStrands;
			buildResultCount( work->config, work->primary, pipeline.bestOffsets,
				pipeline.bestStrands, initial.count );
			initial.consensus = buildConsensus( work->config, work->primary, initial.count );
			refineMotifWithCache( work->config, work->primary, work->control,
				work->primaryCache, work->controlCache, *work->table,
				initial, &(*work->results)[pipelineIdx] );
		} catch ( ... ) {
			(*work->errors)[pipelineIdx] = std::current_exception();
		}
	}
}

}


RefinementSelection selectRefinementThreshold( const std::vector<RefinementSite> &primarySites,
		const std::vector<RefinementSite> &controlSites ) {
	if ( primarySites.empty() || controlSites.empty() ) return RefinementSelection{};
	validatePopulationSizes( primarySites.size(), controlSites.size() );
	const FisherTable table( primarySites.size() + controlSites.size() );
	return selectRefinementThresholdWithTable( primarySites, controlSites, table );
}


void refineMotif( const Config *config, const Dataset *primary, const Dataset *control,
		const RefinementBackground *background, const OutputMotif &startingMotif,
		RefinementResult *result ) {
	validatePopulationSizes( primary->sequences.size(), control->sequences.size() );
	const FisherTable table( primary->sequences.size() + control->sequences.size() );
	RefinementBackgroundCache primaryCache;
	RefinementBackgroundCache controlCache;
	buildRefinementBackgroundCache( config, primary, background, &primaryCache );
	buildRefinementBackgroundCache( config, control, background, &controlCache );
	refineMotifWithCache( config, primary, control, &primaryCache, &controlCache,
		table, startingMotif, result );
}


void refinePipelineMotifs( const Config *config, const Dataset *primary, const Dataset *control,
		const std::vector<PipelineResult> &pipelineResults,
		RefinementBackground *background, std::vector<RefinementResult> &results ) {
	buildRefinementBackground( control, background );
	results.clear();
	results.resize( pipelineResults.size() );
	if ( pipelineResults.empty() ) return;

	// Build once per invocation, outside all candidate and refinement loops.
	validatePopulationSizes( primary->sequences.size(), control->sequences.size() );
	const FisherTable table( primary->sequences.size() + control->sequences.size() );
	RefinementBackgroundCache primaryCache;
	RefinementBackgroundCache controlCache;
	buildRefinementBackgroundCache( config, primary, background, &primaryCache );
	buildRefinementBackgroundCache( config, control, background, &controlCache );

	int threadNum = config->threadNum;
	if ( threadNum == 0 ) threadNum = NUMPIPELINE;
	threadNum = std::max( 1, std::min( threadNum, NUMPIPELINE ) );
	size_t workerNum = std::min( (size_t)threadNum, pipelineResults.size() );
	std::vector<std::exception_ptr> errors( pipelineResults.size() );
	RefinementWork work{};
	work.config = config;
	work.primary = primary;
	work.control = control;
	work.pipelineResults = &pipelineResults;
	work.primaryCache = &primaryCache;
	work.controlCache = &controlCache;
	work.table = &table;
	work.results = &results;
	work.errors = &errors;
	if ( workerNum == 1 ) {
		refinePipelineWorker( &work );
	} else {
		std::vector<std::thread> workers;
		workers.reserve( workerNum );
		try {
			for ( size_t workerIdx = 0; workerIdx < workerNum; workerIdx ++ ) {
				workers.emplace_back( refinePipelineWorker, &work );
			}
		} catch ( ... ) {
			for ( size_t workerIdx = 0; workerIdx < workers.size(); workerIdx ++ ) {
				workers[workerIdx].join();
			}
			throw;
		}
		for ( size_t workerIdx = 0; workerIdx < workers.size(); workerIdx ++ ) {
			workers[workerIdx].join();
		}
	}
	// Propagate errors only after every worker has released the shared tables.
	for ( size_t pipelineIdx = 0; pipelineIdx < errors.size(); pipelineIdx ++ ) {
		if ( errors[pipelineIdx] ) std::rethrow_exception( errors[pipelineIdx] );
	}
}
