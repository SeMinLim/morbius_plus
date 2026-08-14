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
	// Select the best pipeline and store the result
	//--------------------------------------------------------------------------------------------
	printf( "[STEP 4] Selecting and storing the best result is started!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
	int bestPipelineIdx = selectBestPipeline(pipelineResults);
	const vector<uint32_t> &bestOffsets = pipelineResults[bestPipelineIdx].bestOffsets;
	vector<uint32_t> resultCount;
	buildResultCount(&config, &dataset, bestOffsets, resultCount);
	string consensus = buildConsensus(&config, &dataset, resultCount);
	createOutputDirectory(config.outputPrefix);
	writeMotifFASTA(&config, &dataset, bestOffsets);
	writeOffsets(&config, &dataset, bestOffsets);
	writePWM(&config, &dataset, resultCount);
	writeMEME(&config, &dataset, resultCount);
	writeSummary(&config,
		     &dataset,
		     &seedModel,
		     pipelineResults,
		     bestPipelineIdx,
		     consensus,
		     seedElapsedTime + processElapsedTime);
	printf( "[STEP 4] Selecting and storing the best result is done!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	printResult(&config,
		    &dataset,
		    &seedModel,
		    pipelineResults,
		    bestPipelineIdx,
		    consensus,
		    seedElapsedTime + processElapsedTime);

	return 0;
}
