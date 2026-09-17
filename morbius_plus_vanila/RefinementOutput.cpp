#include "Refinement.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

// Compare exact smoothed fractions without overflow in cross products.
static bool sameRefinedPWM( const OutputMotif &first, const OutputMotif &second ) {
	if ( first.count.size() != second.count.size() ) return false;
	uint64_t firstDenominator = (uint64_t)first.siteNum + 4 * PSEUDOCOUNT;
	uint64_t secondDenominator = (uint64_t)second.siteNum + 4 * PSEUDOCOUNT;
	for ( size_t idx = 0; idx < first.count.size(); idx ++ ) {
		uint64_t firstNumerator = (uint64_t)first.count[idx] + PSEUDOCOUNT;
		uint64_t secondNumerator = (uint64_t)second.count[idx] + PSEUDOCOUNT;
		uint64_t firstDivisor = gcd(firstNumerator, firstDenominator);
		uint64_t secondDivisor = gcd(secondNumerator, secondDenominator);
		if ( firstNumerator / firstDivisor != secondNumerator / secondDivisor ||
		     firstDenominator / firstDivisor != secondDenominator / secondDivisor ) return false;
	}
	return true;
}

// Refinement does not alter the original Gibbs score ordering.
void selectRefinedOutputMotifs( const Config *config,
		const vector<PipelineResult> &pipelineResults,
		const vector<RefinementResult> &results, vector<OutputMotif> &motifs ) {
	motifs.clear();
	vector<size_t> order(results.size());
	iota(order.begin(), order.end(), 0);
	for ( const RefinementResult &result : results ) {
		if ( result.motif.pipelineIdx < 0 ||
		     (size_t)result.motif.pipelineIdx >= pipelineResults.size() ) {
			printf( "Invalid pipeline index in refinement result.\n" );
			exit(1);
		}
	}
	sort(order.begin(), order.end(), [&](size_t first, size_t second) {
		int firstPipeline = results[first].motif.pipelineIdx;
		int secondPipeline = results[second].motif.pipelineIdx;
		uint64_t firstScore = pipelineResults[firstPipeline].bestScore;
		uint64_t secondScore = pipelineResults[secondPipeline].bestScore;
		return firstScore != secondScore ? firstScore > secondScore : firstPipeline < secondPipeline;
	});
	for ( size_t resultIdx : order ) {
		if ( (uint64_t)motifs.size() >= config->outputMotifNum ) break;
		const RefinementResult &result = results[resultIdx];
		if ( result.selection.valid == false || result.motif.siteNum == 0 ) continue;
		bool duplicate = false;
		for ( const OutputMotif &motif : motifs ) {
			if ( sameRefinedPWM(result.motif, motif) ) {
				duplicate = true;
				break;
			}
		}
		if ( duplicate == false ) motifs.push_back(result.motif);
	}
}

static ofstream openRefinementFile( const string &filename ) {
	ofstream outputFile(filename);
	if ( outputFile.is_open() == false ) {
		printf( "Unable to create output file: %s\n", filename.c_str() );
		exit(1);
	}
	outputFile << setprecision(17);
	return outputFile;
}

static void closeRefinementFile( ofstream &outputFile, const string &filename ) {
	outputFile.close();
	if ( outputFile.fail() ) {
		printf( "Unable to write output file: %s\n", filename.c_str() );
		exit(1);
	}
}

static void writeRefinementMEMEHeader( ofstream &outputFile ) {
	outputFile << "MEME version 4\n\nALPHABET= ACGT\n\nstrands: +\n\n";
	// Match the existing result MEME so default Tomtom scoring does not depend on
	// whether candidates were exported before or after the final output limit.
	// Actual refinement likelihoods use the separate Control Markov tables.
	outputFile << "Background letter frequencies\n";
	for ( size_t symbol = 0; symbol < 4; symbol ++ ) {
		if ( symbol != 0 ) outputFile << " ";
		outputFile << "ACGT"[symbol] << " " << 0.25;
	}
	outputFile << "\n\n";
}

static void writeRefinementMEMEMotif( ofstream &outputFile, const Config *config,
		const OutputMotif &motif, const vector<double> *pwm, size_t siteNum ) {
	outputFile << "MOTIF MorbiusPlus_pipeline_" << motif.pipelineIdx << "\n";
	outputFile << "letter-probability matrix: alength= 4 w= " << config->motifLength
		   << " nsites= " << siteNum << " E= 0\n";
	double denominator = (double)siteNum + 4 * PSEUDOCOUNT;
	for ( size_t column = 0; column < config->motifLength; column ++ ) {
		for ( size_t symbol = 0; symbol < 4; symbol ++ ) {
			if ( symbol != 0 ) outputFile << " ";
			size_t idx = symbol * config->motifLength + column;
			double probability = pwm != NULL ? (*pwm)[idx] :
				((double)motif.count[idx] + PSEUDOCOUNT) / denominator;
			outputFile << probability;
		}
		outputFile << "\n";
	}
	outputFile << "\n";
}

static string contextString( size_t code, size_t order ) {
	if ( order == 0 ) return "NA";
	string context(order, 'A');
	for ( size_t idx = order; idx > 0; idx -- ) {
		context[idx - 1] = "ACGT"[code % 4];
		code /= 4;
	}
	return context;
}

// Preserve every candidate, including candidates omitted by deduplication or the output limit.
void writeRefinement( const Config *config, const Dataset *primary, const Dataset *control,
		const RefinementBackground *background, const vector<RefinementResult> &results,
		const vector<OutputMotif> &motifs, double elapsedTime ) {
	string controlSource = config->controlFilename.empty() ?
		"Generated third-order Markov DNA in memory from Primary" : config->controlFilename;
	string filename = config->outputPrefix + ".refinement.tsv";
	ofstream outputFile = openRefinementFile(filename);
	outputFile << "Pipeline\tOutputRank\tFitSource\tPrimarySiteNum"
		   << "\tEvaluationPrimaryHits\tEvaluationControlHits\tEvaluationScoreThreshold\tEvaluationLogPvalue"
		   << "\tFitPrimaryHits\tFitControlHits\tFitScoreThreshold\tFitSelectionLogPvalue"
		   << "\tInitialLogPvalue\tAcceptedIterations\tAttemptedIterations\tTerminationReason\n";
	for ( const RefinementResult &result : results ) {
		const OutputMotif &motif = result.motif;
		outputFile << motif.pipelineIdx << "\t";
		size_t rank = 0;
		for ( size_t idx = 0; idx < motifs.size(); idx ++ ) {
			if ( motifs[idx].pipelineIdx == motif.pipelineIdx ) rank = idx + 1;
		}
		if ( rank == 0 ) outputFile << "unreported";
		else outputFile << rank;
		bool supported = result.selection.valid && motif.siteNum > 0;
		outputFile << "\t" << (supported ?
			(result.acceptedIterations == 0 ? "original_gibbs" : "refinement") : "unsupported")
			<< "\t" << motif.siteNum << "\t";
		if ( result.selection.valid ) {
			outputFile << result.selection.primarySiteNum << "\t" << result.selection.controlSiteNum
				<< "\t" << result.selection.scoreThreshold << "\t" << result.selection.logPvalue;
		} else outputFile << "NA\tNA\tNA\tNA";
		outputFile << "\t";
		if ( result.fittingSelection.valid && result.acceptedIterations > 0 ) {
			outputFile << result.fittingSelection.primarySiteNum << "\t" << result.fittingSelection.controlSiteNum
				<< "\t" << result.fittingSelection.scoreThreshold << "\t" << result.fittingSelection.logPvalue;
		} else outputFile << "NA\tNA\tNA\tNA";
		outputFile << "\t" << result.initialLogPvalue << "\t" << result.acceptedIterations
			   << "\t" << result.iterations << "\t" << result.terminationReason << "\n";
	}
	closeRefinementFile(outputFile, filename);

	filename = config->outputPrefix + ".refinement_sites.tsv";
	outputFile = openRefinementFile(filename);
	outputFile << "Pipeline\tSequenceIdx\tSequenceName\tSitePresent\tOffset\tStrand\n";
	for ( const RefinementResult &result : results ) {
		const OutputMotif &motif = result.motif;
		for ( size_t seqIdx = 0; seqIdx < primary->sequences.size(); seqIdx ++ ) {
			bool present = seqIdx < motif.sitePresent.size() && motif.sitePresent[seqIdx] != 0;
			outputFile << motif.pipelineIdx << "\t" << seqIdx << "\t" << primary->names[seqIdx]
				   << "\t" << (present ? 1 : 0) << "\t";
			if ( present ) {
				outputFile << motif.offsets[seqIdx] << "\t"
					   << (motif.strands[seqIdx] == STRAND_REVERSE ? '-' : '+') << "\n";
			} else outputFile << "NA\tNA\n";
		}
	}
	closeRefinementFile(outputFile, filename);

	filename = config->outputPrefix + ".refinement_scoring.meme";
	outputFile = openRefinementFile(filename);
	writeRefinementMEMEHeader(outputFile);
	outputFile << "# These preceding PWMs selected the fitted sites for the corresponding candidates.\n";
	outputFile << "# nsites describes the sites originally used to construct each preceding PWM.\n";
	outputFile << "# The fitted candidate PWM and its evaluated enrichment are stored separately.\n\n";
	for ( const RefinementResult &result : results ) {
		if ( result.selection.valid && result.motif.siteNum > 0 && result.scoringPWM.empty() == false ) {
			writeRefinementMEMEMotif(outputFile, config, result.motif, &result.scoringPWM, result.scoringSiteNum);
		}
	}
	closeRefinementFile(outputFile, filename);

	filename = config->outputPrefix + ".refinement_candidates.meme";
	outputFile = openRefinementFile(filename);
	writeRefinementMEMEHeader(outputFile);
	for ( const RefinementResult &result : results ) {
		if ( result.selection.valid && result.motif.siteNum > 0 ) {
			writeRefinementMEMEMotif(outputFile, config, result.motif, NULL, result.motif.siteNum);
		}
	}
	closeRefinementFile(outputFile, filename);

	filename = config->outputPrefix + ".refinement_background.tsv";
	outputFile = openRefinementFile(filename);
	outputFile << "# Control=" << controlSource << "\n"
		   << "# DNA Markov order=" << REFINEMENTBACKGROUNDORDER << "; both Control strands; prior 1/4^k per length-k tuple, then conditional normalization\n"
		   << "Order\tContext\tBase\tLog2ConditionalProbability\tConditionalProbability\n";
	for ( size_t order = 0; order < background->logProbability.size(); order ++ ) {
		for ( size_t idx = 0; idx < background->logProbability[order].size(); idx ++ ) {
			double logProbability = background->logProbability[order][idx];
			outputFile << order << "\t" << contextString(idx / 4, order) << "\t" << "ACGT"[idx % 4]
				   << "\t" << logProbability << "\t" << exp2(logProbability) << "\n";
		}
	}
	closeRefinementFile(outputFile, filename);

	filename = config->outputPrefix + ".refinement_summary.txt";
	outputFile = openRefinementFile(filename);
	outputFile << "Method: STREME-inspired post-Gibbs DNA ZOOPS refinement; not an exact STREME implementation.\n"
		   << "Primary: " << config->inputFilename << "\nControl: " << controlSource << "\n"
		   << "Primary sequences: " << primary->sequences.size() << "\nControl sequences: " << control->sequences.size() << "\n"
		   << "Primary sequence length: " << primary->sequenceLength << "\nControl sequence length: " << control->sequenceLength << "\n"
		   << "Motif width: " << config->motifLength << "\nCandidate pipelines: " << results.size() << "\n"
		   << "Reported motifs: " << motifs.size() << "\nMaximum refinement iterations: " << REFINEMENTMAXITERATIONS << "\n"
		   << "Gibbs agreement threshold: " << config->scoreThreshold << " (unchanged; separate from site-score thresholds)\n"
		   << "Background: Control only, both strands, Markov order " << REFINEMENTBACKGROUNDORDER << ", computed once and reused.\n"
		   << "Background smoothing: prior 1/4^k per length-k tuple (total prior mass 1 per tuple length), then conditional normalization.\n"
		   << "Background context: site-local, starting at order 0 for the first oriented base.\n"
		   << "Scan: best log2(PWM/background) site per original sequence over forward strand only.\n"
		   << "Threshold: complete score-tie groups, inclusive score >= threshold; minimum one-sided Fisher p-value.\n"
		   << "PWM fit: at most one passing Primary site per sequence; Control sites never enter PWM counts.\n"
		   << "Pseudocount: " << PSEUDOCOUNT << " per base; PWM probability=(count+1)/(PrimarySiteNum+4).\n"
		   << "Acceptance: re-evaluate the reconstructed candidate PWM; retain it only if enrichment strictly improves.\n"
		   << "Evaluation fields: optimal threshold, hit counts, and natural-log p-value measured on the final candidate PWM.\n"
		   << "Fit fields and scoring MEME: preceding PWM threshold and sites used to reconstruct that candidate.\n"
		   << "Fit sites need not equal the final candidate's evaluation hits; both stages are recorded separately.\n"
		   << "Original Gibbs fitting sites are retained when no refinement step improves an enriched starting PWM.\n"
		   << "Unsupported candidates have no reported PWM; all candidates remain in refinement diagnostics.\n"
		   << "Output order: unchanged original Gibbs bestScore descending, pipeline index ascending on ties.\n"
		   << "Deduplication: exact smoothed PWM probabilities after refinement; no reverse-complement or similarity clustering.\n"
		   << "Refinement elapsed seconds: " << elapsedTime << "\n"
		   << "Program elapsed includes input, seed initialization, Gibbs, and refinement; it stops before motif selection and output. External timing covers the whole process.\n"
		   << "No reference motifs or Tomtom scores are used in refinement.\n";
	if ( config->controlFilename.empty() ) {
		outputFile << "Control generation: third-order Markov model from all forward Primary 1-mer through 4-mer observations; no pseudocounts on observed rows.\n"
			   << "Control backoff: unobserved transition rows use shorter suffixes; an empty order-zero row uses uniform probabilities.\n"
			   << "Control sampling: orders 0, 1 and 2 for the first three bases, then order 3; floor(2^32 * cumulative / total) integer boundaries.\n"
			   << "Control RNG: per-sequence xoshiro128+, initialized with SplitMix64(state = fixed seed + zero-based sequence index); fixed seed "
			   << DEFAULTCONTROLSEED << "; independent of Gibbs seed and thread scheduling.\n"
			   << "Control generated once in memory; same sequence counts and lengths as Primary; counting, model preparation and generation included in program elapsed.\n";
	}
	outputFile << "All MEME exports retain uniform letter frequencies for the existing Tomtom evaluation convention; refinement uses the separate Control Markov background.\n";
	closeRefinementFile(outputFile, filename);
}
