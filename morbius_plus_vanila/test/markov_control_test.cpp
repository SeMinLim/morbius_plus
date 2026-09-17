#include "MorbiusPlus.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

static void require( bool condition, const string &message ) {
	if ( condition ) return;
	cerr << "Markov Control test failed: " << message << "\n";
	exit(1);
}

// Enumerate substrings independently of the streaming packed-context helper.
static ControlCounts manualCounts( const vector<string> &sequences ) {
	ControlCounts counts = {};
	const size_t starts[] = {0, 4, 20, 84};
	for ( const string &sequence : sequences ) {
		for ( size_t length = 1; length <= 4; length ++ ) {
			for ( size_t position = 0; position + length <= sequence.size(); position ++ ) {
				size_t code = 0;
				bool valid = true;
				for ( size_t idx = 0; idx < length; idx ++ ) {
					size_t symbol = string("ACGT").find(sequence[position + idx]);
					if ( symbol == string::npos ) { valid = false; break; }
					code = code * 4 + symbol;
				}
				if ( valid ) counts.kmer[starts[length - 1] + code] ++;
			}
		}
	}
	return counts;
}

static void sameCounts( const ControlCounts &actual, const ControlCounts &expected ) {
	for ( size_t idx = 0; idx < 340; idx ++ ) {
		require(actual.kmer[idx] == expected.kmer[idx], "k-mer count differs at flat index " + to_string(idx));
	}
}

static void checkCounting( const string &filename ) {
	const vector<string> sequences = {"ACGTAC", "T", "AAAANCGT", "GCXTA"};
	ControlCounts actual = {};
	for ( const string &sequence : sequences ) {
		uint32_t context = 0;
		size_t length = 0;
		for ( char base : sequence ) {
			size_t symbol = string("ACGT").find(base);
			countControlBase(&actual, &context, &length, symbol == string::npos ? -1 : (int)symbol);
		}
	}
	sameCounts(actual, manualCounts(sequences));
	uint32_t context = 7;
	size_t length = 3;
	countControlBase(&actual, &context, &length, 4);
	require(context == 0 && length == 0, "invalid high symbol did not reset history");

	ofstream output(filename, ios::binary);
	output << ">first\r\naCgT\r\ntAcG\r\n>second\r\nTT\r\n aa cc GG\r\n";
	output.close();
	require(!output.fail(), "cannot create wrapped FASTA fixture");
	Dataset parsed;
	configureAlphabet(ALPHABET_DNA, &parsed);
	ControlCounts parsedCounts = {};
	readFASTA(filename, &parsed, &parsedCounts);
	require(parsed.sequences == vector<string>({"ACGTTACG", "TTAACCGG"}), "FASTA normalization changed");
	sameCounts(parsedCounts, manualCounts(parsed.sequences));
	Dataset legacy;
	configureAlphabet(ALPHABET_DNA, &legacy);
	readFASTA(filename, &legacy);
	require(legacy.sequences == parsed.sequences && legacy.names == parsed.names,
		"optional counting changed FASTA input");
}

static void rowEquals( const ControlTransitions &table, size_t row,
		uint64_t first, uint64_t second, uint64_t third ) {
	require(table.cumulative[row][0] == first && table.cumulative[row][1] == second &&
		table.cumulative[row][2] == third, "unexpected CDF at row " + to_string(row));
}

static void checkTransitions( void ) {
	ControlCounts counts = {};
	ControlTransitions table;
	buildControlTransitions(&counts, &table);
	for ( size_t row = 0; row < 85; row ++ ) rowEquals(table, row, 1073741824ULL, 2147483648ULL, 3221225472ULL);
	counts = manualCounts({"AAAAC", "CAAAG"});
	buildControlTransitions(&counts, &table);
	// AAA and CAA share the same suffix AA but have different observed outgoing bases.
	rowEquals(table, 21, 1431655765ULL, 2863311530ULL, 4294967296ULL);
	rowEquals(table, 21 + 16, 4294967296ULL, 4294967296ULL, 4294967296ULL);
	for ( size_t edge = 0; edge < 3; edge ++ ) {
		require(table.cumulative[21 + 32][edge] == table.cumulative[5][edge],
			"unseen GAA did not back off to observed AA");
		require(table.cumulative[21 + 7][edge] == table.cumulative[0][edge],
			"unseen ACT did not recursively back off CT and T to order zero");
	}
	counts = {};
	counts.kmer[1] = 1;
	buildControlTransitions(&counts, &table);
	rowEquals(table, 0, 0, 4294967296ULL, 4294967296ULL);
	counts = {};
	counts.kmer[0] = 18446744073709551613ULL;
	counts.kmer[1] = counts.kmer[2] = 1;
	buildControlTransitions(&counts, &table);
	rowEquals(table, 0, 4294967295ULL, 4294967295ULL, 4294967296ULL);
	counts = {};
	counts.kmer[0] = 9223372036854775808ULL;
	counts.kmer[3] = 1;
	buildControlTransitions(&counts, &table);
	rowEquals(table, 0, 4294967295ULL, 4294967295ULL, 4294967295ULL);
}

static void checkGeneration( void ) {
	Config config = {};
	config.alphabetMode = ALPHABET_DNA;
	config.randomSeed = 1;
	config.threadNum = 1;
	Dataset primary;
	configureAlphabet(ALPHABET_DNA, &primary);
	const size_t lengths[] = {1, 2, 3, 31, 32, 33, 100};
	const string pattern = "ACGTTGCATGTCGCATGATGCATGAGAGCT";
	for ( size_t idx = 0; idx < 7; idx ++ ) {
		primary.names.push_back("sequence_" + to_string(idx));
		string sequence;
		while ( sequence.size() < lengths[idx] ) sequence += pattern;
		primary.sequences.push_back(sequence.substr(0, lengths[idx]));
	}
	primary.sequenceLength = lengths[0];
	ControlCounts counts = manualCounts(primary.sequences);
	const vector<string> expected = {
		"G", "AT", "CTA", "AGAGAGAGAGAGAGCTACGTTGCATGAGAGC",
		"ACGTTGCATGATGTCGCATGAGCTACGTTGCA", "CTACGTTGCATGATGCATGTCGCATGTCGCATG",
		"GTTGCATGTCGCATGTCGCATGTCGCATGCATGAGAGCTACGTTGCATGAGAGCTACGTTGCATGTCGCATGCATGTCGCATGAGCTACGTTGCATGTCG"
	};
	require(DEFAULTCONTROLSEED == 1, "fixed Control seed changed");
	Dataset control;
	generateMarkovControl(&config, &primary, &counts, &control);
	require(control.sequences == expected, "fixed order-three sampling vector changed");
	require(control.names == primary.names && control.sequenceLength == primary.sequenceLength &&
		control.alphabet == "ACGT" && control.alphabetSize == 4, "Control shape or metadata changed");
	for ( size_t idx = 0; idx < 7; idx ++ ) {
		require(control.sequences[idx].size() == lengths[idx], "individual sequence length changed");
		for ( char base : control.sequences[idx] ) require(control.symbolMap[(unsigned char)base] == (int)string("ACGT").find(base),
			"DNA symbol map is incorrect");
	}
	for ( int threads : {4, 0, 16, 64} ) {
		config.threadNum = threads;
		config.randomSeed = 741;
		Dataset repeated;
		generateMarkovControl(&config, &primary, &counts, &repeated);
		require(repeated.sequences == expected, "Control depends on thread count, Gibbs seed or retained RNG state");
	}
	RandomGenerator withControl, withoutControl;
	initializeRandomGenerator(&withControl, 741);
	initializeRandomGenerator(&withoutControl, 741);
	generateMarkovControl(&config, &primary, &counts, &control);
	for ( size_t idx = 0; idx < 100; idx ++ ) require(randomWord(&withControl) == randomWord(&withoutControl),
		"Control generation changed Gibbs RNG state");
	counts = {};
	counts.kmer[2] = 1;
	generateMarkovControl(&config, &primary, &counts, &control);
	for ( size_t idx = 0; idx < 7; idx ++ ) require(control.sequences[idx] == string(lengths[idx], 'G'),
		"zero-probability bases were sampled or unseen contexts failed suffix backoff");
}

int main( int argc, char **argv ) {
	require(argc == 2, "expected temporary FASTA path");
	checkCounting(argv[1]);
	checkTransitions();
	checkGeneration();
	cout << "Validated streaming counts, true order-three transitions, exact CDF boundaries, parser normalization, and independent parallel sampling.\n";
	return 0;
}
