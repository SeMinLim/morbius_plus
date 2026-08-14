#include "SyntheticDatasetGenerator.h"

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <random>
#include <string>
#include <vector>
using namespace std;


// Validate all binding-site sequences
size_t validateBindingSites( vector<FastaRecord> &bindingSites ) {
	size_t bindingSiteLength = bindingSites[0].sequence.size();
	if ( bindingSiteLength < 3 || bindingSiteLength > DNA_SEQUENCE_LENGTH ) {
		printf( "Unsupported binding-site length: %lu\n", (unsigned long)bindingSiteLength );
		fflush( stdout );
		exit(1);
	}

	for ( size_t siteIdx = 0; siteIdx < bindingSites.size(); siteIdx ++ ) {
		if ( bindingSites[siteIdx].sequence.size() != bindingSiteLength ) {
			printf( "Binding sites must have an identical length. Record %lu has length %lu instead of %lu.\n",
				(unsigned long)siteIdx,
				(unsigned long)bindingSites[siteIdx].sequence.size(),
				(unsigned long)bindingSiteLength );
			fflush( stdout );
			exit(1);
		}

		for ( size_t position = 0; position < bindingSiteLength; position ++ ) {
			bindingSites[siteIdx].sequence[position] =
				toUpperSymbol(bindingSites[siteIdx].sequence[position]);
			if ( bindingSites[siteIdx].sequence[position] == 'U' ) {
				bindingSites[siteIdx].sequence[position] = 'T';
			}
			if ( isValidDNA(bindingSites[siteIdx].sequence[position]) == false ) {
				printf( "Unsupported symbol in binding-site record %lu at position %lu.\n",
					(unsigned long)siteIdx,
					(unsigned long)position );
				fflush( stdout );
				exit(1);
			}
		}
	}

	return bindingSiteLength;
}

// Chop input background records into non-overlapping 1000-base sequences
void buildBackgroundPool( const vector<FastaRecord> &inputRecords,
			  vector<string> &backgroundPool ) {
	for ( size_t recordIdx = 0; recordIdx < inputRecords.size(); recordIdx ++ ) {
		const string &sequence = inputRecords[recordIdx].sequence;
		size_t chunkNum = sequence.size() / DNA_SEQUENCE_LENGTH;
		for ( size_t chunkIdx = 0; chunkIdx < chunkNum; chunkIdx ++ ) {
			backgroundPool.push_back(sequence.substr(chunkIdx * DNA_SEQUENCE_LENGTH,
								   DNA_SEQUENCE_LENGTH));
		}
	}

	if ( backgroundPool.empty() ) {
		printf( "The background FASTA does not contain a complete 1000-base sequence.\n" );
		fflush( stdout );
		exit(1);
	}
}

// Return one amino acid for a DNA codon
// TAA and TAG map to Q, and TGA maps to W, to keep the 20-symbol protein alphabet
char translateCodon( char base0, char base1, char base2 ) {
	if ( base0 == 'T' ) {
		if ( base1 == 'T' ) {
			if ( base2 == 'T' || base2 == 'C' ) return 'F';
			return 'L';
		}
		if ( base1 == 'C' ) return 'S';
		if ( base1 == 'A' ) {
			if ( base2 == 'T' || base2 == 'C' ) return 'Y';
			return 'Q';
		}
		if ( base2 == 'T' || base2 == 'C' ) return 'C';
		return 'W';
	}

	if ( base0 == 'C' ) {
		if ( base1 == 'T' ) return 'L';
		if ( base1 == 'C' ) return 'P';
		if ( base1 == 'A' ) {
			if ( base2 == 'T' || base2 == 'C' ) return 'H';
			return 'Q';
		}
		return 'R';
	}

	if ( base0 == 'A' ) {
		if ( base1 == 'T' ) {
			if ( base2 == 'G' ) return 'M';
			return 'I';
		}
		if ( base1 == 'C' ) return 'T';
		if ( base1 == 'A' ) {
			if ( base2 == 'T' || base2 == 'C' ) return 'N';
			return 'K';
		}
		if ( base2 == 'T' || base2 == 'C' ) return 'S';
		return 'R';
	}

	if ( base1 == 'T' ) return 'V';
	if ( base1 == 'C' ) return 'A';
	if ( base1 == 'A' ) {
		if ( base2 == 'T' || base2 == 'C' ) return 'D';
		return 'E';
	}
	return 'G';
}

// Translate complete codons from one binding site
string translateBindingSite( const string &bindingSite ) {
	size_t translatedLength = bindingSite.size() / 3;
	string proteinMotif;
	proteinMotif.reserve(translatedLength);

	for ( size_t position = 0; position + 2 < bindingSite.size(); position += 3 ) {
		char aminoAcid = translateCodon(bindingSite[position],
						bindingSite[position + 1],
						bindingSite[position + 2]);
		proteinMotif.push_back(aminoAcid);
	}

	return proteinMotif;
}

// Build translated protein motifs from the DNA binding-site source
void buildProteinMotifs( const vector<FastaRecord> &bindingSites,
			 vector<string> &proteinMotifs ) {
	proteinMotifs.resize(bindingSites.size());
	for ( size_t siteIdx = 0; siteIdx < bindingSites.size(); siteIdx ++ ) {
		proteinMotifs[siteIdx] = translateBindingSite(bindingSites[siteIdx].sequence);
		if ( proteinMotifs[siteIdx].empty() ||
		     proteinMotifs[siteIdx].size() > PROTEIN_SEQUENCE_LENGTH ) {
			printf( "Unsupported translated motif length for binding-site record %lu.\n",
				(unsigned long)siteIdx );
			fflush( stdout );
			exit(1);
		}
	}
}

// Build dataset specifications
void buildDatasetSpecs( bool testMode, DatasetSpec specs[DATASET_NUM] ) {
	if ( testMode ) {
		specs[0] = {1, 8, 2};
		specs[1] = {2, 16, 4};
		specs[2] = {3, 32, 8};
	} else {
		specs[0] = {1, 32768, 1024};
		specs[1] = {2, 65536, 2048};
		specs[2] = {3, 131072, 4096};
	}
}

// Select unique source indices
void selectUniqueIndices( uint64_t sourceNum,
			  uint64_t selectedNum,
			  mt19937_64 &randomGenerator,
			  vector<uint32_t> &selectedIndices ) {
	if ( selectedNum > sourceNum ) {
		printf( "Requested %lu unique entries from a source containing only %lu entries.\n",
			(unsigned long)selectedNum,
			(unsigned long)sourceNum );
		fflush( stdout );
		exit(1);
	}

	vector<uint32_t> sourceIndices(sourceNum);
	for ( uint64_t i = 0; i < sourceNum; i ++ ) sourceIndices[i] = (uint32_t)i;
	shuffle(sourceIndices.begin(), sourceIndices.end(), randomGenerator);
	selectedIndices.assign(sourceIndices.begin(), sourceIndices.begin() + selectedNum);
}
