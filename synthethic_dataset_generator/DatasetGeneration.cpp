#include "SyntheticDatasetGenerator.h"

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <random>
#include <string>
#include <vector>
using namespace std;


// Write one FASTA record
void writeFASTARecord( ofstream &output, uint64_t sequenceID, const string &sequence ) {
	output << ">" << sequenceID << "\n";
	if ( FASTA_LINE_WIDTH == 0 ) {
		output << sequence << "\n";
	} else {
		for ( size_t position = 0; position < sequence.size(); position += FASTA_LINE_WIDTH ) {
			output << sequence.substr(position, FASTA_LINE_WIDTH) << "\n";
		}
	}
}

// Generate one DNA dataset and its codon-translated protein dataset
void generateDatasetPair( const Config *config,
			  const DatasetSpec *spec,
			  const vector<string> &backgroundPool,
			  const vector<FastaRecord> &bindingSites ) {
	uint64_t datasetSeed = config->randomSeed ^
			       ((uint64_t)spec->datasetID * 0x9e3779b97f4a7c15ULL);
	mt19937_64 randomGenerator(datasetSeed);

	vector<uint32_t> backgroundIndices;
	vector<uint32_t> motifPoolIndices;
	selectUniqueIndices(backgroundPool.size(),
			    spec->sequenceNum,
			    randomGenerator,
			    backgroundIndices);
	selectUniqueIndices(bindingSites.size(),
			    spec->motifPoolNum,
			    randomGenerator,
			    motifPoolIndices);

	string datasetID = to_string(spec->datasetID);
	string dnaFilename = config->outputDirectory + "/DATASET_DNA_" + datasetID + ".fasta";
	string dnaGroundTruthFilename = config->outputDirectory +
					"/DATASET_DNA_" + datasetID + "_ground_truth.tsv";
	string proteinFilename = config->outputDirectory +
				 "/DATASET_PROTEIN_" + datasetID + ".fasta";
	string proteinGroundTruthFilename = config->outputDirectory +
					    "/DATASET_PROTEIN_" + datasetID + "_ground_truth.tsv";

	ofstream dnaOutput(dnaFilename.c_str());
	ofstream dnaGroundTruthOutput(dnaGroundTruthFilename.c_str());
	ofstream proteinOutput(proteinFilename.c_str());
	ofstream proteinGroundTruthOutput(proteinGroundTruthFilename.c_str());
	if ( dnaOutput.is_open() == false ||
	     dnaGroundTruthOutput.is_open() == false ||
	     proteinOutput.is_open() == false ||
	     proteinGroundTruthOutput.is_open() == false ) {
		printf( "Unable to create dataset output files for dataset %d.\n", spec->datasetID );
		fflush( stdout );
		exit(1);
	}

	dnaGroundTruthOutput << "sequence_id\tbackground_id\tbinding_site_id\timplanted_offset\timplanted_length\timplanted_sequence\n";
	proteinGroundTruthOutput << "sequence_id\tbinding_site_id\tsource_dna_offset\tsource_dna_length\ttranslated_offset\ttranslated_length\ttranslated_sequence\n";

	uniform_int_distribution<uint64_t> motifPoolDistribution(0, spec->motifPoolNum - 1);
	uint64_t progressInterval = spec->sequenceNum / 32;
	if ( progressInterval < 1 ) progressInterval = 1;

	for ( uint64_t sequenceID = 0; sequenceID < spec->sequenceNum; sequenceID ++ ) {
		uint32_t backgroundID = backgroundIndices[sequenceID];
		uint32_t bindingSiteID = motifPoolIndices[motifPoolDistribution(randomGenerator)];
		const string &bindingSite = bindingSites[bindingSiteID].sequence;

		// Generate the DNA sequence and record the exact implanted position
		string dnaSequence = backgroundPool[backgroundID];
		normalizeDNA(dnaSequence, randomGenerator);
		uniform_int_distribution<uint64_t> dnaOffsetDistribution(
			0,
			DNA_SEQUENCE_LENGTH - bindingSite.size()
		);
		uint64_t dnaOffset = dnaOffsetDistribution(randomGenerator);
		dnaSequence.replace(dnaOffset, bindingSite.size(), bindingSite);
		if ( dnaSequence.substr(dnaOffset, bindingSite.size()) != bindingSite ) {
			printf( "DNA implantation validation failed for dataset %d, sequence %lu.\n",
				spec->datasetID,
				(unsigned long)sequenceID );
			fflush( stdout );
			exit(1);
		}

		writeFASTARecord(dnaOutput, sequenceID, dnaSequence);
		dnaGroundTruthOutput << sequenceID << "\t"
					 << backgroundID << "\t"
					 << bindingSiteID << "\t"
					 << dnaOffset << "\t"
					 << bindingSite.size() << "\t"
					 << bindingSite << "\n";

		// Translate the generated DNA sequence with the standard DNA codon table
		string proteinSequence;
		string translatedSequence;
		int64_t translatedOffset = -1;
		translateDNASequence(dnaSequence,
				     dnaOffset,
				     bindingSite.size(),
				     proteinSequence,
				     &translatedOffset,
				     translatedSequence);
		writeFASTARecord(proteinOutput, sequenceID, proteinSequence);
		proteinGroundTruthOutput << sequenceID << "\t"
					     << bindingSiteID << "\t"
					     << dnaOffset << "\t"
					     << bindingSite.size() << "\t"
					     << translatedOffset << "\t"
					     << translatedSequence.size() << "\t"
					     << translatedSequence << "\n";

		if ( (sequenceID + 1) % progressInterval == 0 ||
		     sequenceID + 1 == spec->sequenceNum ) {
			printf( "[STEP 3] Dataset %d generation is processing...[%lu/%lu]\n",
				spec->datasetID,
				(unsigned long)(sequenceID + 1),
				(unsigned long)spec->sequenceNum );
			fflush( stdout );
		}
	}

	dnaOutput.close();
	dnaGroundTruthOutput.close();
	proteinOutput.close();
	proteinGroundTruthOutput.close();

	printf( "[STEP 3] Dataset %d generation is done!\n", spec->datasetID );
	printf( "DNA Dataset Size       : %lu x %d\n",
		(unsigned long)spec->sequenceNum,
		DNA_SEQUENCE_LENGTH );
	printf( "Protein Dataset Size   : %lu x %d\n",
		(unsigned long)spec->sequenceNum,
		PROTEIN_SEQUENCE_LENGTH );
	printf( "Motif Pool Size        : %lu\n", (unsigned long)spec->motifPoolNum );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}
