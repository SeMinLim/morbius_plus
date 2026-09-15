#ifndef REFINEMENT_H
#define REFINEMENT_H

#include "MorbiusPlus.h"

#define REFINEMENTMAXITERATIONS 20
#define REFINEMENTBACKGROUNDORDER 2

typedef struct RefinementBackground {
	// Conditional log2 probabilities, indexed by order and (context * 4 + base).
	std::vector<std::vector<double>> logProbability;
	std::vector<double> letterFrequency;
}RefinementBackground;

typedef struct RefinementSite {
	double score;
	uint32_t offset;
	uint8_t strand;
}RefinementSite;

typedef struct RefinementSelection {
	bool valid = false;
	double scoreThreshold = 0.0;
	double logPvalue = 0.0; // Natural logarithm of the one-sided Fisher p-value.
	size_t primarySiteNum = 0;
	size_t controlSiteNum = 0;
}RefinementSelection;

typedef struct RefinementResult {
	OutputMotif motif;
	RefinementSelection selection; // Re-evaluated on the final output PWM.
	RefinementSelection fittingSelection; // Selected sites used to construct that PWM.
	double initialLogPvalue = 0.0;
	size_t iterations = 0;
	size_t acceptedIterations = 0;
	std::string terminationReason;
	// The fitting threshold describes sites selected using this preceding PWM.
	// Empty when the original Gibbs fit is retained without an accepted refinement.
	std::vector<double> scoringPWM;
	size_t scoringSiteNum = 0;
	std::vector<double> acceptedLogPvalues;
}RefinementResult;

void buildRefinementBackground( const Dataset *control, RefinementBackground *background );
double calculateRefinementBackgroundLogProbability( const Config *config,
		const Dataset *dataset, const RefinementBackground *background,
		const std::string &sequence, uint32_t offset, uint8_t strand );
void scanRefinementSites( const Config *config, const Dataset *dataset,
		const RefinementBackground *background, const std::vector<double> &pwm,
		std::vector<RefinementSite> &sites );
double fisherLogPvalue( size_t primaryHits, size_t controlHits,
		size_t primaryNum, size_t controlNum );
RefinementSelection selectRefinementThreshold( const std::vector<RefinementSite> &primarySites,
		const std::vector<RefinementSite> &controlSites );
void refineMotif( const Config *config, const Dataset *primary, const Dataset *control,
		const RefinementBackground *background, const OutputMotif &startingMotif,
		RefinementResult *result );
void refinePipelineMotifs( const Config *config, const Dataset *primary, const Dataset *control,
		const std::vector<PipelineResult> &pipelineResults,
		RefinementBackground *background, std::vector<RefinementResult> &results );
void selectRefinedOutputMotifs( const Config *config,
		const std::vector<PipelineResult> &pipelineResults,
		const std::vector<RefinementResult> &results, std::vector<OutputMotif> &motifs );
void writeRefinement( const Config *config, const Dataset *primary, const Dataset *control,
		const RefinementBackground *background, const std::vector<RefinementResult> &results,
		const std::vector<OutputMotif> &motifs, double elapsedTime );

#endif
