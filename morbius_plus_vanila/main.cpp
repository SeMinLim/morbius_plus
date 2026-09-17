#include "MorbiusPlus.h"
#include "Refinement.h"

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>
using namespace std;


// Main
int main( int argc, char **argv ) {
	double programStartTime = timeChecker();
	Config config;
	parseArguments(argc, argv, &config);

	Dataset dataset;
	configureAlphabet(config.alphabetMode, &dataset);

	//--------------------------------------------------------------------------------------------
	// [STEP 1]
	// Read and validate the input sequence dataset
	//--------------------------------------------------------------------------------------------
	printf( "---------------------------------------------------------------------\n" );
	printf( "[STEP 1] Reading sequence FASTA file is started!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
	bool generateControl = config.alphabetMode == ALPHABET_DNA && config.controlFilename.empty();
	ControlCounts controlCounts;
	readFASTA(config.inputFilename, &dataset, generateControl ? &controlCounts : NULL);
	validateWorkload(&config, &dataset);
	Dataset control;
	bool refinementEnabled = config.alphabetMode == ALPHABET_DNA;
	if ( refinementEnabled ) {
		if ( config.controlFilename.empty() ) {
			generateMarkovControl(&config, &dataset, &controlCounts, &control);
		} else {
			configureAlphabet(config.alphabetMode, &control);
			readFASTA(config.controlFilename, &control);
		}
		validateWorkload(&config, &control);
		if ( control.sequenceLength != dataset.sequenceLength ) {
			fprintf(stderr, "Fisher refinement requires the same sequence length in Primary and Control.\n");
			return 1;
		}
	}
	if ( config.maxUpdateNum == 0 ) {
		config.maxUpdateNum = (uint64_t)DEFAULTMAXSWEEPNUM * (uint64_t)dataset.sequences.size();
	}
	printf( "[STEP 1] Reading sequence FASTA file is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	//--------------------------------------------------------------------------------------------
	// [STEP 2]
	// Build the fixed-rate support-guided initialization model
	//--------------------------------------------------------------------------------------------
	printf( "[STEP 2] Building the support-guided initialization model is started!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
	double seedStartTime = timeChecker();
	vector<SeedModel> seedModels;
	buildSeedModels(&config, &dataset, seedModels);
	double seedElapsedTime = timeChecker() - seedStartTime;
	printf( "[STEP 2] Building the support-guided initialization model is done!\n" );
	printf( "Seed Initialization Time: %.8f\n", seedElapsedTime );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	//--------------------------------------------------------------------------------------------
	// [STEP 3]
	// Run independent Morbius+ Gibbs pipelines
	//--------------------------------------------------------------------------------------------
	printf( "[STEP 3] Running Morbius+ Gibbs pipelines is started!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
	double processStartTime = timeChecker();
	vector<PipelineResult> pipelineResults;
	runPipelines(&config, &dataset, seedModels, pipelineResults);
	double processElapsedTime = timeChecker() - processStartTime;
	printf( "[STEP 3] Running Morbius+ Gibbs pipelines is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	// Separate postprocessing: all Gibbs results and their RNG state remain unchanged.
	RefinementBackground refinementBackground;
	vector<RefinementResult> refinementResults;
	double refinementElapsedTime = 0.0;
	if ( refinementEnabled ) {
		printf( "[Refinement] Refining all %lu Gibbs candidates with Primary and Control.\n",
			(unsigned long)pipelineResults.size() );
		fflush(stdout);
		double refinementStartTime = timeChecker();
		refinePipelineMotifs(&config, &dataset, &control, pipelineResults,
			&refinementBackground, refinementResults);
		refinementElapsedTime = timeChecker() - refinementStartTime;
	}
	// Stop before motif selection and all result-file output.
	// DNA includes input reading and any Control generation; protein retains seed + Gibbs timing.
	double elapsedTime = refinementEnabled ? timeChecker() - programStartTime :
		seedElapsedTime + processElapsedTime;

	//--------------------------------------------------------------------------------------------
	// [STEP 4]
	// Select unique output motifs and store the results
	//--------------------------------------------------------------------------------------------
	if ( config.outputMotifNum == 1 ) {
		printf( "[STEP 4] Selecting and storing the best result is started!\n" );
	} else {
		printf( "[STEP 4] Selecting and storing output motifs is started!\n" );
	}
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
	vector<OutputMotif> outputMotifs;
	if ( refinementEnabled ) {
		selectRefinedOutputMotifs(&config, pipelineResults, refinementResults, outputMotifs);
	} else {
		selectOutputMotifs(&config, &dataset, pipelineResults, outputMotifs);
	}
	createOutputDirectory(config.outputPrefix);
	writeInitialization(&config, &dataset, seedModels, pipelineResults);
	writeMotifFASTA(&config, &dataset, pipelineResults, outputMotifs);
	writeOffsets(&config, &dataset, pipelineResults, outputMotifs);
	writePWM(&config, &dataset, outputMotifs);
	writeMEME(&config, &dataset, outputMotifs);
	if ( refinementEnabled ) {
		writeRefinement(&config, &dataset, &control, &refinementBackground,
			refinementResults, outputMotifs, refinementElapsedTime);
	}
	writeSummary(&config,
		     &dataset,
		     seedModels,
		     pipelineResults,
		     outputMotifs,
		     elapsedTime);
	if ( config.outputMotifNum == 1 ) {
		printf( "[STEP 4] Selecting and storing the best result is done!\n" );
	} else {
		printf( "[STEP 4] Selecting and storing output motifs is done!\n" );
	}
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	printResult(&config,
		    &dataset,
		    seedModels,
		    pipelineResults,
		    outputMotifs,
		    elapsedTime);

	return 0;
}
