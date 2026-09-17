#include "MorbiusPlus.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

static void require( bool condition, const string &message ) {
	if ( condition ) return;
	cerr << "Uniform Control test failed: " << message << "\n";
	exit(1);
}

int main( void ) {
	Dataset primary;
	configureAlphabet(ALPHABET_DNA, &primary);
	const size_t lengths[] = {31, 32, 33, 100};
	primary.sequenceLength = lengths[0];
	for ( size_t idx = 0; idx < 4; idx ++ ) {
		primary.names.push_back("primary_" + to_string(idx));
		primary.sequences.push_back(string(lengths[idx], 'G'));
	}
	const string expected =
		"CAATATCCGAAACGAGATGTCTGAGGAACACGTCGCATGTGTAGCCGCCAGGCTAGTGGTGTTGG"
		"TCCCCCCGATATGTTGTGTGAGGTACGAGTTTGAACGATGAACGTGTAACGGCAGCAATCATCCG"
		"TGCCTGCAAACACTAGCTACCCTGTGCATCAAAGGAAACCCAAACGTTTTTGAACTACTAATCCGG";
	require(DEFAULTCONTROLSEED == 1, "documented fixed seed changed");
	Dataset control;
	generateUniformControl(&primary, &control);
	require(control.names == primary.names && control.sequences.size() == primary.sequences.size(),
		"sequence names or count do not match Primary");
	require(control.sequenceLength == primary.sequenceLength && control.alphabet == "ACGT" &&
		control.alphabetSize == 4, "generated Dataset metadata is incorrect");
	size_t position = 0;
	vector<size_t> counts(4, 0);
	for ( size_t idx = 0; idx < 4; idx ++ ) {
		require(control.sequences[idx].size() == lengths[idx], "individual sequence length changed");
		require(control.sequences[idx] == expected.substr(position, lengths[idx]),
			"fixed vector or 32-base word reuse across sequence boundaries changed");
		for ( char base : control.sequences[idx] ) {
			size_t symbol = control.alphabet.find(base);
			require(symbol < 4 && control.symbolMap[(unsigned char)base] == (int)symbol,
				"invalid generated symbol or DNA mapping");
			counts[symbol] ++;
		}
		position += lengths[idx];
	}
	require(position == expected.size() && counts == vector<size_t>({53, 48, 50, 45}),
		"fixed vector base counts changed");
	Dataset repeated;
	generateUniformControl(&primary, &repeated);
	require(repeated.sequences == control.sequences, "generation retained hidden state between calls");
	for ( string &sequence : primary.sequences ) sequence.assign(sequence.size(), 'A');
	generateUniformControl(&primary, &repeated);
	require(repeated.sequences == control.sequences, "uniform generator depends on Primary base frequencies");
	// A live Gibbs stream must be identical whether Control generation occurs or not.
	RandomGenerator withControl, withoutControl;
	initializeRandomGenerator(&withControl, 741);
	initializeRandomGenerator(&withoutControl, 741);
	for ( size_t idx = 0; idx < 100; idx ++ ) {
		generateUniformControl(&primary, &repeated);
		require(randomWord(&withControl) == randomWord(&withoutControl), "Control generation changed Gibbs RNG state");
	}
	cout << "Validated fixed Control vector, per-sequence shape, word boundaries, and independent RNG.\n";
	return 0;
}
