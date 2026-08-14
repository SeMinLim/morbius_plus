#include "SyntheticDatasetGenerator.h"

#include <sys/stat.h>
#include <sys/time.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <random>
#include <string>
#include <vector>
using namespace std;


// Time checker
double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)(tv.tv_sec) + (double)(tv.tv_usec) / 1000000.0;
}

// Print command usage
void printUsage( const char *programName ) {
	printf( "Usage: %s --background <upstream5000.fa> --binding-sites <MA0007.2_sites.fasta> --output <directory> [--seed <value>] [--test-mode]\n",
		programName );
	printf( "\n" );
	printf( "Required inputs\n" );
	printf( "  --background      FASTA file containing the upstream background sequences\n" );
	printf( "  --binding-sites   FASTA file containing equal-length DNA binding-site sequences\n" );
	printf( "  --output          Output directory for DNA1-3, Protein1-3, and ground truth\n" );
	printf( "\n" );
	printf( "Optional inputs\n" );
	printf( "  --seed            Deterministic random seed (default: 1)\n" );
	printf( "  --test-mode       Generate 8, 16, and 32 sequences instead of production sizes\n" );
	fflush( stdout );
}

// Parse command-line arguments
void parseArguments( int argc, char **argv, Config *config ) {
	config->randomSeed = 1;
	config->testMode = false;

	for ( int argIdx = 1; argIdx < argc; argIdx ++ ) {
		string argument = argv[argIdx];
		if ( argument == "--background" ) {
			if ( argIdx + 1 >= argc ) {
				printf( "Missing value after --background.\n" );
				exit(1);
			}
			config->backgroundFilename = argv[++ argIdx];
		} else if ( argument == "--binding-sites" ) {
			if ( argIdx + 1 >= argc ) {
				printf( "Missing value after --binding-sites.\n" );
				exit(1);
			}
			config->bindingSiteFilename = argv[++ argIdx];
		} else if ( argument == "--output" ) {
			if ( argIdx + 1 >= argc ) {
				printf( "Missing value after --output.\n" );
				exit(1);
			}
			config->outputDirectory = argv[++ argIdx];
		} else if ( argument == "--seed" ) {
			if ( argIdx + 1 >= argc ) {
				printf( "Missing value after --seed.\n" );
				exit(1);
			}
			config->randomSeed = strtoull(argv[++ argIdx], NULL, 10);
		} else if ( argument == "--test-mode" ) {
			config->testMode = true;
		} else if ( argument == "--help" || argument == "-h" ) {
			printUsage(argv[0]);
			exit(0);
		} else {
			printf( "Unknown argument: %s\n", argument.c_str() );
			printUsage(argv[0]);
			exit(1);
		}
	}

	if ( config->backgroundFilename.empty() ||
	     config->bindingSiteFilename.empty() ||
	     config->outputDirectory.empty() ) {
		printUsage(argv[0]);
		exit(1);
	}
}

// Create an output directory
void createDirectory( const string &directory ) {
	if ( mkdir(directory.c_str(), 0755) != 0 && errno != EEXIST ) {
		printf( "Unable to create output directory: %s\n", directory.c_str() );
		fflush( stdout );
		exit(1);
	}
}

// Convert a character to upper case without locale dependency
char toUpperSymbol( char symbol ) {
	if ( symbol >= 'a' && symbol <= 'z' ) return (char)(symbol - 'a' + 'A');
	return symbol;
}

// Read a FASTA file
void readFASTA( const string &filename, vector<FastaRecord> &records ) {
	ifstream input(filename.c_str());
	if ( input.is_open() == false ) {
		printf( "File not found: %s\n", filename.c_str() );
		fflush( stdout );
		exit(1);
	}

	string line;
	FastaRecord record;
	while ( getline(input, line) ) {
		if ( line.empty() ) continue;
		if ( line[line.size() - 1] == '\r' ) line.pop_back();
		if ( line.empty() ) continue;

		if ( line[0] == '>' ) {
			if ( record.sequence.empty() == false ) {
				records.push_back(record);
				record = FastaRecord();
			}
			record.name = line.substr(1);
		} else {
			for ( size_t i = 0; i < line.size(); i ++ ) {
				char symbol = toUpperSymbol(line[i]);
				if ( symbol != ' ' && symbol != '\t' ) record.sequence.push_back(symbol);
			}
		}
	}

	if ( record.sequence.empty() == false ) records.push_back(record);
	input.close();

	if ( records.empty() ) {
		printf( "No FASTA sequences were found in: %s\n", filename.c_str() );
		fflush( stdout );
		exit(1);
	}
}

// Check whether a DNA symbol is valid
bool isValidDNA( char symbol ) {
	return symbol == 'A' || symbol == 'C' || symbol == 'G' || symbol == 'T';
}

// Replace unsupported DNA symbols with deterministic random bases
void normalizeDNA( string &sequence, mt19937_64 &randomGenerator ) {
	uniform_int_distribution<int> baseDistribution(0, 3);
	static const char DNA_BASES[5] = "ACGT";

	for ( size_t i = 0; i < sequence.size(); i ++ ) {
		sequence[i] = toUpperSymbol(sequence[i]);
		if ( sequence[i] == 'U' ) sequence[i] = 'T';
		if ( isValidDNA(sequence[i]) == false ) {
			sequence[i] = DNA_BASES[baseDistribution(randomGenerator)];
		}
	}
}
