#include "MorbiusPlus.h"

#include <stdio.h>

#include <string>
#include <vector>
using namespace std;


// Main
int main( int argc, char **argv ) {
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
	readFASTA(config.inputFilename, &dataset);
	validateWorkload(&config, &dataset);
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
	SeedModel seedModel;
	buildSeedModel(&config, &dataset, &seedModel);
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
	runPipelines(&config, &dataset, &seedModel, pipelineResults);
	double processElapsedTime = timeChecker() - processStartTime;
	printf( "[STEP 3] Running Morbius+ Gibbs pipelines is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

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
	selectOutputMotifs(&config, &dataset, pipelineResults, outputMotifs);
	createOutputDirectory(config.outputPrefix);
	writeMotifFASTA(&config, &dataset, pipelineResults, outputMotifs);
	writeOffsets(&config, &dataset, pipelineResults, outputMotifs);
	writePWM(&config, &dataset, outputMotifs);
	writeMEME(&config, &dataset, outputMotifs);
	writeSummary(&config,
		     &dataset,
		     &seedModel,
		     pipelineResults,
		     outputMotifs,
		     seedElapsedTime + processElapsedTime);
	if ( config.outputMotifNum == 1 ) {
		printf( "[STEP 4] Selecting and storing the best result is done!\n" );
	} else {
		printf( "[STEP 4] Selecting and storing output motifs is done!\n" );
	}
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	printResult(&config,
		    &dataset,
		    &seedModel,
		    pipelineResults,
		    outputMotifs,
		    seedElapsedTime + processElapsedTime);

	return 0;
}
