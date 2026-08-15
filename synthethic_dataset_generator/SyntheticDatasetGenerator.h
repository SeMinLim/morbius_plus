#ifndef SYNTHETIC_DATASET_GENERATOR_H
#define SYNTHETIC_DATASET_GENERATOR_H

#include <stdint.h>

#include <fstream>
#include <random>
#include <string>
#include <vector>


#define DNA_SEQUENCE_LENGTH 1000
#define PROTEIN_SEQUENCE_LENGTH 300
#define DATASET_NUM 3
#define FASTA_LINE_WIDTH 0


typedef struct Config {
	std::string backgroundFilename;
	std::string bindingSiteFilename;
	std::string outputDirectory;
	uint64_t randomSeed;
	bool testMode;
}Config;

typedef struct DatasetSpec {
	int datasetID;
	uint64_t sequenceNum;
	uint64_t motifPoolNum;
}DatasetSpec;

typedef struct FastaRecord {
	std::string name;
	std::string sequence;
}FastaRecord;


double timeChecker( void );
void printUsage( const char *programName );
void parseArguments( int argc, char **argv, Config *config );
void createDirectory( const std::string &directory );
char toUpperSymbol( char symbol );
void readFASTA( const std::string &filename, std::vector<FastaRecord> &records );
bool isValidDNA( char symbol );
void normalizeDNA( std::string &sequence, std::mt19937_64 &randomGenerator );

size_t validateBindingSites( std::vector<FastaRecord> &bindingSites );
void buildBackgroundPool( const std::vector<FastaRecord> &inputRecords,
			  std::vector<std::string> &backgroundPool );
char translateCodon( char base0, char base1, char base2 );
void translateDNASequence( const std::string &dnaSequence,
			   uint64_t implantedOffset,
			   uint64_t implantedLength,
			   std::string &proteinSequence,
			   int64_t *translatedOffset,
			   std::string &translatedSequence );
void buildDatasetSpecs( bool testMode, DatasetSpec specs[DATASET_NUM] );
void selectUniqueIndices( uint64_t sourceNum,
			  uint64_t selectedNum,
			  std::mt19937_64 &randomGenerator,
			  std::vector<uint32_t> &selectedIndices );

void writeFASTARecord( std::ofstream &output,
		       uint64_t sequenceID,
		       const std::string &sequence );
void generateDatasetPair( const Config *config,
			  const DatasetSpec *spec,
			  const std::vector<std::string> &backgroundPool,
			  const std::vector<FastaRecord> &bindingSites );

#endif
