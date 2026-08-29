#include "MorbiusPlus.h"

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <string>
using namespace std;


// Elapsed time checker
double timeChecker( void ) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

// Integer power
uint64_t integerPower( uint64_t base, size_t exponent ) {
	uint64_t result = 1;
	for ( size_t i = 0; i < exponent; i ++ ) result = result * base;
	return result;
}

// Rotate a 32-bit value left
static inline uint32_t rotateLeft32( uint32_t value, int shift ) {
	return (value << shift) | (value >> (32 - shift));
}

// SplitMix64
uint64_t splitMix64( uint64_t *state ) {
	uint64_t result = (*state += 0x9e3779b97f4a7c15ULL);
	result = (result ^ (result >> 30)) * 0xbf58476d1ce4e5b9ULL;
	result = (result ^ (result >> 27)) * 0x94d049bb133111ebULL;
	return result ^ (result >> 31);
}

// Initialize xoshiro128+
void initializeRandomGenerator( RandomGenerator *randomGenerator, uint64_t seed ) {
	uint64_t seedState = seed;
	for ( int i = 0; i < 4; i ++ ) {
		uint64_t value = splitMix64(&seedState);
		randomGenerator->state[i] = (uint32_t)(value ^ (value >> 32));
	}

	if ( randomGenerator->state[0] == 0 && randomGenerator->state[1] == 0 &&
	     randomGenerator->state[2] == 0 && randomGenerator->state[3] == 0 ) {
		randomGenerator->state[0] = 1;
	}
}

// Generate a 32-bit random word with xoshiro128+
uint32_t randomWord( RandomGenerator *randomGenerator ) {
	uint32_t result = randomGenerator->state[0] + randomGenerator->state[3];
	uint32_t shiftValue = randomGenerator->state[1] << 9;

	randomGenerator->state[2] = randomGenerator->state[2] ^ randomGenerator->state[0];
	randomGenerator->state[3] = randomGenerator->state[3] ^ randomGenerator->state[1];
	randomGenerator->state[1] = randomGenerator->state[1] ^ randomGenerator->state[2];
	randomGenerator->state[0] = randomGenerator->state[0] ^ randomGenerator->state[3];
	randomGenerator->state[2] = randomGenerator->state[2] ^ shiftValue;
	randomGenerator->state[3] = rotateLeft32(randomGenerator->state[3], 11);

	return result;
}

// Generate a 24-bit random value
uint32_t randomQ24( RandomGenerator *randomGenerator ) {
	return randomWord(randomGenerator) >> (32 - RANDOMFRACTIONBITS);
}

// Generate a uniform floating-point random value
double randomUnit( RandomGenerator *randomGenerator ) {
	return (double)randomQ24(randomGenerator) / (double)(1U << RANDOMFRACTIONBITS);
}

// Generate a bounded integer random value
uint32_t randomBounded( RandomGenerator *randomGenerator, uint32_t bound ) {
	if ( bound == 0 ) return 0;
	return randomWord(randomGenerator) % bound;
}

// Print usage
void printUsage( const char *programName ) {
	printf( "Usage: %s --input <FASTA> --output <PREFIX> --alphabet <dna|protein> --motif-length <N> [Options]\n", programName );
	printf( "\n" );
	printf( "Options:\n" );
	printf( "  --max-updates <N>      Maximum updates per pipeline [default: %d x sequence count]\n", DEFAULTMAXSWEEPNUM );
	printf( "  --score-threshold <F>  Normalized agreement threshold in [0, 1] [default: %.2f]\n", DEFAULTSCORETHRESHOLD );
	printf( "  --seed <N>             Random seed [default: %d]\n", DEFAULTSEED );
	printf( "  --threads <N>          Concurrent CPU threads [default: %d]\n", NUMPIPELINE );
	printf( "  --help                 Print this message\n" );
}

// Parse an unsigned integer
uint64_t parseUnsignedInteger( const char *value, const char *name ) {
	if ( value[0] == '-' ) {
		printf( "Invalid value for %s: %s\n", name, value );
		exit(1);
	}
	char *end = NULL;
	unsigned long long result = strtoull(value, &end, 10);
	if ( end == value || *end != '\0' ) {
		printf( "Invalid value for %s: %s\n", name, value );
		exit(1);
	}
	return (uint64_t)result;
}

// Parse a floating-point value
double parseFloatingPoint( const char *value, const char *name ) {
	char *end = NULL;
	double result = strtod(value, &end);
	if ( end == value || *end != '\0' ) {
		printf( "Invalid value for %s: %s\n", name, value );
		exit(1);
	}
	return result;
}

// Parse command-line arguments
void parseArguments( int argc, char **argv, Config *config ) {
	config->alphabetMode = -1;
	config->motifLength = 0;
	config->maxUpdateNum = 0;
	config->scoreThreshold = DEFAULTSCORETHRESHOLD;
	config->randomSeed = DEFAULTSEED;
	config->threadNum = 0;

	for ( int i = 1; i < argc; i ++ ) {
		if ( string(argv[i]) == "--input" && i + 1 < argc ) {
			config->inputFilename = argv[++i];
		} else if ( string(argv[i]) == "--output" && i + 1 < argc ) {
			config->outputPrefix = argv[++i];
		} else if ( string(argv[i]) == "--alphabet" && i + 1 < argc ) {
			string alphabetName = argv[++i];
			if ( alphabetName == "dna" || alphabetName == "DNA" ) config->alphabetMode = ALPHABET_DNA;
			else if ( alphabetName == "protein" || alphabetName == "PROTEIN" ) config->alphabetMode = ALPHABET_PROTEIN;
			else {
				printf( "Unsupported alphabet: %s\n", alphabetName.c_str() );
				exit(1);
			}
		} else if ( string(argv[i]) == "--motif-length" && i + 1 < argc ) {
			config->motifLength = (size_t)parseUnsignedInteger(argv[++i], "--motif-length");
		} else if ( string(argv[i]) == "--max-updates" && i + 1 < argc ) {
			config->maxUpdateNum = parseUnsignedInteger(argv[++i], "--max-updates");
		} else if ( string(argv[i]) == "--score-threshold" && i + 1 < argc ) {
			config->scoreThreshold = parseFloatingPoint(argv[++i], "--score-threshold");
		} else if ( string(argv[i]) == "--seed" && i + 1 < argc ) {
			config->randomSeed = parseUnsignedInteger(argv[++i], "--seed");
		} else if ( string(argv[i]) == "--threads" && i + 1 < argc ) {
			config->threadNum = (int)parseUnsignedInteger(argv[++i], "--threads");
		} else if ( string(argv[i]) == "--help" ) {
			printUsage(argv[0]);
			exit(0);
		} else {
			printf( "Unknown or incomplete argument: %s\n", argv[i] );
			printUsage(argv[0]);
			exit(1);
		}
	}

	if ( config->inputFilename.empty() || config->outputPrefix.empty() ||
	     config->alphabetMode < 0 || config->motifLength == 0 ) {
		printUsage(argv[0]);
		exit(1);
	}
	if ( config->scoreThreshold < 0.0 || config->scoreThreshold > 1.0 ) {
		printf( "The score threshold must be in [0, 1].\n" );
		exit(1);
	}
	if ( config->threadNum < 0 ) {
		printf( "The number of threads cannot be negative.\n" );
		exit(1);
	}
}

// Configure the sequence alphabet
void configureAlphabet( int alphabetMode, Dataset *dataset ) {
	dataset->symbolMap.assign(256, -1);
	if ( alphabetMode == ALPHABET_DNA ) {
		dataset->alphabet = "ACGT";
	} else {
		dataset->alphabet = "ACDEFGHIKLMNPQRSTVWY";
	}
	dataset->alphabetSize = (int)dataset->alphabet.size();

	for ( int i = 0; i < dataset->alphabetSize; i ++ ) {
		unsigned char symbol = (unsigned char)dataset->alphabet[i];
		dataset->symbolMap[symbol] = i;
	}
}

// Store a FASTA sequence
void storeSequence( const string &name, const string &sequence, Dataset *dataset ) {
	if ( sequence.empty() ) return;

	string sequenceUpper = sequence;
	for ( size_t i = 0; i < sequenceUpper.size(); i ++ ) {
		unsigned char symbol = (unsigned char)sequenceUpper[i];
		if ( symbol >= 'a' && symbol <= 'z' ) sequenceUpper[i] = (char)(symbol - 'a' + 'A');
		if ( dataset->symbolMap[(unsigned char)sequenceUpper[i]] < 0 ) {
			printf( "Unsupported symbol '%c' in sequence %s at position %lu.\n",
				sequenceUpper[i],
				name.c_str(),
				(unsigned long)i
			);
			exit(1);
		}
	}

	if ( dataset->sequences.empty() ) dataset->sequenceLength = sequenceUpper.size();
	else if ( sequenceUpper.size() != dataset->sequenceLength ) {
		printf( "All sequences must have the same length. Sequence %s has length %lu instead of %lu.\n",
			name.c_str(),
			(unsigned long)sequenceUpper.size(),
			(unsigned long)dataset->sequenceLength
		);
		exit(1);
	}

	dataset->names.push_back(name.empty() ? to_string(dataset->names.size()) : name);
	dataset->sequences.push_back(sequenceUpper);
}

// Read a FASTA dataset
void readFASTA( const string &filename, Dataset *dataset ) {
	ifstream inputFile(filename);
	if ( inputFile.is_open() == false ) {
		printf( "File not found: %s\n", filename.c_str() );
		exit(1);
	}

	string line;
	string sequenceName;
	string sequence;
	while ( getline(inputFile, line) ) {
		if ( line.empty() ) continue;
		if ( line[0] == '>' ) {
			storeSequence(sequenceName, sequence, dataset);
			sequenceName = line.substr(1);
			sequence.clear();
		} else {
			for ( size_t i = 0; i < line.size(); i ++ ) {
				unsigned char symbol = (unsigned char)line[i];
				if ( symbol != ' ' && symbol != '\t' && symbol != '\r' ) sequence.push_back((char)symbol);
			}
		}
	}
	storeSequence(sequenceName, sequence, dataset);
	inputFile.close();

	if ( dataset->sequences.empty() ) {
		printf( "No sequence was found in %s.\n", filename.c_str() );
		exit(1);
	}
}

// Validate the workload
void validateWorkload( const Config *config, const Dataset *dataset ) {
	if ( config->motifLength > dataset->sequenceLength ) {
		printf( "The motif length cannot exceed the sequence length.\n" );
		exit(1);
	}
	if ( config->alphabetMode == ALPHABET_DNA && config->motifLength < 4 ) {
		printf( "DNA motif length must be at least 4.\n" );
		exit(1);
	}
	if ( config->alphabetMode == ALPHABET_PROTEIN && config->motifLength < 3 ) {
		printf( "Protein motif length must be at least 3.\n" );
		exit(1);
	}
}
