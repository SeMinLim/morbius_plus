#include "Refinement.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>


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


double fisherLogPvalue( size_t primaryHits, size_t controlHits,
		size_t primaryNum, size_t controlNum ) {
	if ( primaryHits > primaryNum || controlHits > controlNum ) {
		throw std::invalid_argument( "Fisher hit counts exceed their sequence populations." );
	}
	validatePopulationSizes( primaryNum, controlNum );
	FisherTable table( primaryNum + controlNum );
	return (double)evaluateFisherLogPvalue( table, primaryHits, controlHits, primaryNum, controlNum );
}


RefinementSelection selectRefinementThreshold( const std::vector<RefinementSite> &primarySites,
		const std::vector<RefinementSite> &controlSites ) {
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
	FisherTable table( total );
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


void refineMotif( const Config *config, const Dataset *primary, const Dataset *control,
		const RefinementBackground *background, const OutputMotif &startingMotif,
		RefinementResult *result ) {
	OutputMotif initial = startingMotif;
	*result = RefinementResult{};
	result->motif = initial;
	result->motif.refined = true;
	std::vector<double> pwm;
	rebuildRefinementPWM( config, result->motif, pwm );
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


void refinePipelineMotifs( const Config *config, const Dataset *primary, const Dataset *control,
		const std::vector<PipelineResult> &pipelineResults,
		RefinementBackground *background, std::vector<RefinementResult> &results ) {
	buildRefinementBackground( control, background );
	results.clear();
	results.resize( pipelineResults.size() );
	for ( size_t pipelineIdx = 0; pipelineIdx < pipelineResults.size(); pipelineIdx++ ) {
		OutputMotif initial{};
		initial.pipelineIdx = (int)pipelineIdx;
		initial.refined = true;
		initial.siteNum = primary->sequences.size();
		initial.sitePresent.assign( primary->sequences.size(), 1 );
		initial.offsets = pipelineResults[pipelineIdx].bestOffsets;
		initial.strands = pipelineResults[pipelineIdx].bestStrands;
		buildResultCount( config, primary, pipelineResults[pipelineIdx].bestOffsets,
			pipelineResults[pipelineIdx].bestStrands, initial.count );
		initial.consensus = buildConsensus( config, primary, initial.count );
		refineMotif( config, primary, control, background, initial, &results[pipelineIdx] );
	}
}
