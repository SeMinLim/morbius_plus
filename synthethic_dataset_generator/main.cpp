#include "SyntheticDatasetGenerator.h"

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>
using namespace std;


// Main
int main( int argc, char **argv ) {
	Config config;
	parseArguments(argc, argv, &config);
	createDirectory(config.outputDirectory);

	//--------------------------------------------------------------------------------------------
	// [STEP 1]
	// Read source FASTA files
	//--------------------------------------------------------------------------------------------
	printf( "---------------------------------------------------------------------\n" );
	printf( "[STEP 1] Reading source FASTA files started!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	double readStart = timeChecker();
	vector<FastaRecord> backgroundRecords;
	vector<FastaRecord> bindingSites;
	readFASTA(config.backgroundFilename, backgroundRecords);
	readFASTA(config.bindingSiteFilename, bindingSites);
	size_t bindingSiteLength = validateBindingSites(bindingSites);
	double readFinish = timeChecker();

	printf( "[STEP 1] Reading source FASTA files is done!\n" );
	printf( "Background Records     : %lu\n", (unsigned long)backgroundRecords.size() );
	printf( "Binding-Site Records   : %lu\n", (unsigned long)bindingSites.size() );
	printf( "Binding-Site Length    : %lu\n", (unsigned long)bindingSiteLength );
	printf( "Elapsed Time           : %.8f sec\n", readFinish - readStart );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	//--------------------------------------------------------------------------------------------
	// [STEP 2]
	// Prepare the background and protein-motif pools
	//--------------------------------------------------------------------------------------------
	printf( "[STEP 2] Preparing source pools started!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	double poolStart = timeChecker();
	vector<string> backgroundPool;
	vector<string> proteinMotifs;
	buildBackgroundPool(backgroundRecords, backgroundPool);
	buildProteinMotifs(bindingSites, proteinMotifs);
	double poolFinish = timeChecker();

	DatasetSpec specs[DATASET_NUM];
	buildDatasetSpecs(config.testMode, specs);
	uint64_t maximumSequenceNum = specs[DATASET_NUM - 1].sequenceNum;
	uint64_t maximumMotifPoolNum = specs[DATASET_NUM - 1].motifPoolNum;
	if ( backgroundPool.size() < maximumSequenceNum ) {
		printf( "Insufficient background sequences: %lu available, %lu required.\n",
			(unsigned long)backgroundPool.size(),
			(unsigned long)maximumSequenceNum );
		fflush( stdout );
		exit(1);
	}
	if ( bindingSites.size() < maximumMotifPoolNum ) {
		printf( "Insufficient binding sites: %lu available, %lu required.\n",
			(unsigned long)bindingSites.size(),
			(unsigned long)maximumMotifPoolNum );
		fflush( stdout );
		exit(1);
	}

	printf( "[STEP 2] Preparing source pools is done!\n" );
	printf( "1000-Base Backgrounds  : %lu\n", (unsigned long)backgroundPool.size() );
	printf( "Protein Motif Length   : %lu\n", (unsigned long)proteinMotifs[0].size() );
	printf( "Elapsed Time           : %.8f sec\n", poolFinish - poolStart );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	//--------------------------------------------------------------------------------------------
	// [STEP 3]
	// Generate DNA1-3, Protein1-3, and exact ground truth
	//--------------------------------------------------------------------------------------------
	printf( "[STEP 3] Synthetic dataset generation started!\n" );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	double generationStart = timeChecker();
	for ( int datasetIdx = 0; datasetIdx < DATASET_NUM; datasetIdx ++ ) {
		generateDatasetPair(&config,
				    &specs[datasetIdx],
				    backgroundPool,
				    bindingSites,
				    proteinMotifs);
	}
	double generationFinish = timeChecker();

	printf( "[STEP 3] Synthetic dataset generation is done!\n" );
	printf( "Output Directory       : %s\n", config.outputDirectory.c_str() );
	printf( "Random Seed            : %lu\n", (unsigned long)config.randomSeed );
	printf( "Elapsed Time           : %.8f sec\n", generationFinish - generationStart );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );

	return 0;
}
