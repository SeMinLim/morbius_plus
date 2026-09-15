#include "MorbiusPlus.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
using namespace std;

static size_t outputSiteNum( const Dataset *dataset, const OutputMotif &motif ) {
	return motif.refined ? motif.siteNum : dataset->sequences.size();
}

static const vector<uint32_t> &outputOffsets( const vector<PipelineResult> &results,
					     const OutputMotif &motif ) {
	return motif.refined ? motif.offsets : results[motif.pipelineIdx].bestOffsets;
}

static const vector<uint8_t> &outputStrands( const vector<PipelineResult> &results,
					   const OutputMotif &motif ) {
	return motif.refined ? motif.strands : results[motif.pipelineIdx].bestStrands;
}

static void closeOutput( ofstream &outputFile, const string &filename ) {
	outputFile.close();
	if ( outputFile.fail() ) {
		printf( "Unable to write output file: %s\n", filename.c_str() );
		exit(1);
	}
}


// Build result counts without pseudocount
void buildResultCount( const Config *config,
		       const Dataset *dataset,
		       const vector<uint32_t> &offsets,
		       const vector<uint8_t> &strands,
		       vector<uint32_t> &count ) {
	count.assign((size_t)dataset->alphabetSize * config->motifLength, 0);
	for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx ++ ) {
		for ( size_t column = 0; column < config->motifLength; column ++ ) {
			int symbol = getSiteSymbol(config, dataset, dataset->sequences[seqIdx],
						   offsets[seqIdx], strands[seqIdx], column);
			count[(size_t)symbol * config->motifLength + column] ++;
		}
	}
}

// Build the consensus subsequence
string buildConsensus( const Config *config,
		       const Dataset *dataset,
		       const vector<uint32_t> &count ) {
	string consensus(config->motifLength, dataset->alphabet[0]);
	for ( size_t column = 0; column < config->motifLength; column ++ ) {
		uint32_t maximum = 0;
		int maximumSymbol = 0;
		for ( int symbol = 0; symbol < dataset->alphabetSize; symbol ++ ) {
			uint32_t value = count[(size_t)symbol * config->motifLength + column];
			if ( value > maximum ) {
				maximum = value;
				maximumSymbol = symbol;
			}
		}
		consensus[column] = dataset->alphabet[maximumSymbol];
	}
	return consensus;
}

// Select the highest-scoring unique motifs from completed pipelines
void selectOutputMotifs( const Config *config,
			 const Dataset *dataset,
			 const vector<PipelineResult> &pipelineResults,
			 vector<OutputMotif> &outputMotifs ) {
	outputMotifs.clear();
	if ( pipelineResults.empty() ) {
		printf( "No pipeline result is available for output.\n" );
		fflush( stdout );
		exit(1);
	}

	// Preserve the original single-motif Max Filter path
	if ( config->outputMotifNum == 1 ) {
		OutputMotif outputMotif;
		outputMotif.pipelineIdx = selectBestPipeline(pipelineResults);
		buildResultCount(config,
				 dataset,
				 pipelineResults[outputMotif.pipelineIdx].bestOffsets,
				 pipelineResults[outputMotif.pipelineIdx].bestStrands,
				 outputMotif.count);
		outputMotif.consensus = buildConsensus(config, dataset, outputMotif.count);
		outputMotifs.push_back(outputMotif);
		return;
	}

	vector<uint8_t> pipelineProcessed(pipelineResults.size(), 0);
	size_t processedNum = 0;
	while ( processedNum < pipelineResults.size() &&
		(uint64_t)outputMotifs.size() < config->outputMotifNum ) {
		int bestPipelineIdx = -1;
		for ( int pipelineIdx = 0; pipelineIdx < (int)pipelineResults.size(); pipelineIdx ++ ) {
			if ( pipelineProcessed[pipelineIdx] ) continue;
			if ( bestPipelineIdx < 0 ||
			     pipelineResults[pipelineIdx].bestScore >
			     pipelineResults[bestPipelineIdx].bestScore ) {
				bestPipelineIdx = pipelineIdx;
			}
		}

		pipelineProcessed[bestPipelineIdx] = 1;
		processedNum ++;

		OutputMotif outputMotif;
		outputMotif.pipelineIdx = bestPipelineIdx;
		buildResultCount(config,
				 dataset,
				 pipelineResults[bestPipelineIdx].bestOffsets,
				 pipelineResults[bestPipelineIdx].bestStrands,
				 outputMotif.count);

		bool duplicate = false;
		for ( size_t motifIdx = 0; motifIdx < outputMotifs.size(); motifIdx ++ ) {
			if ( outputMotifs[motifIdx].count == outputMotif.count ) {
				duplicate = true;
				break;
			}
		}
		if ( duplicate ) continue;

		outputMotif.consensus = buildConsensus(config, dataset, outputMotif.count);
		outputMotifs.push_back(outputMotif);
	}
}

// Create an output directory
void createOutputDirectory( const string &outputPrefix ) {
	filesystem::path outputPath(outputPrefix);
	filesystem::path parentPath = outputPath.parent_path();
	if ( parentPath.empty() == false ) {
		error_code errorCode;
		filesystem::create_directories(parentPath, errorCode);
		if ( errorCode ) {
			printf( "Unable to create output directory: %s\n", parentPath.string().c_str() );
			exit(1);
		}
	}
}

// Write per-pipeline seed assignments and the offsets captured before Gibbs updates
void writeInitialization( const Config *config,
			  const Dataset *dataset,
			  const vector<SeedModel> &seedModels,
			  const vector<PipelineResult> &pipelineResults ) {
	string seedFilename = config->outputPrefix + ".seeds.tsv";
	ofstream seedFile(seedFilename);
	if ( seedFile.is_open() == false ) {
		printf( "Unable to create output file: %s\n", seedFilename.c_str() );
		exit(1);
	}
	seedFile << setprecision(17);
	seedFile << "Pipeline\tSeedValid\tSeedWord\tSeedCode\tSeedLength\tSeedSupport"
		 << "\tExpectedSeedSupport\tSeedRank\tAnchorOffset\tAnchorRank\n";
	for ( size_t pipelineIdx = 0; pipelineIdx < seedModels.size(); pipelineIdx ++ ) {
		const SeedModel &seedModel = seedModels[pipelineIdx];
		seedFile << pipelineIdx << "\t" << (seedModel.valid ? "Yes" : "No") << "\t";
		if ( seedModel.seedString.empty() ) {
			seedFile << "NA\tNA\t" << seedModel.seedLength << "\t0\tNA\tNA\tNA\tNA\n";
			continue;
		}
		seedFile << seedModel.seedString << "\t" << seedModel.seedCode << "\t"
			 << seedModel.seedLength << "\t" << seedModel.seedSupport << "\t"
			 << seedModel.seedExpectedSupport << "\t" << seedModel.seedRank << "\t";
		if ( seedModel.valid ) {
			seedFile << seedModel.anchorOffset << "\t" << seedModel.anchorRank << "\n";
		} else {
			seedFile << "NA\tNA\n";
		}
	}
	seedFile.close();
	if ( seedFile.fail() ) {
		printf( "Unable to write output file: %s\n", seedFilename.c_str() );
		exit(1);
	}

	string offsetFilename = config->outputPrefix + ".initial_offsets.tsv";
	ofstream offsetFile(offsetFilename);
	if ( offsetFile.is_open() == false ) {
		printf( "Unable to create output file: %s\n", offsetFilename.c_str() );
		exit(1);
	}
	offsetFile << "Sequence\tName";
	for ( size_t pipelineIdx = 0; pipelineIdx < pipelineResults.size(); pipelineIdx ++ ) {
		offsetFile << "\tPipeline_" << pipelineIdx;
	}
	offsetFile << "\n";
	for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx ++ ) {
		offsetFile << seqIdx << "\t" << dataset->names[seqIdx];
		for ( size_t pipelineIdx = 0; pipelineIdx < pipelineResults.size(); pipelineIdx ++ ) {
			offsetFile << "\t" << pipelineResults[pipelineIdx].initialOffsets[seqIdx];
		}
		offsetFile << "\n";
	}
	offsetFile.close();
	if ( offsetFile.fail() ) {
		printf( "Unable to write output file: %s\n", offsetFilename.c_str() );
		exit(1);
	}
}

// Build a site in its selected motif orientation
static string buildSiteString( const Config *config,
			       const Dataset *dataset,
			       const string &sequence,
			       uint32_t offset,
			       uint8_t strand ) {
	string site(config->motifLength, dataset->alphabet[0]);
	for ( size_t column = 0; column < config->motifLength; column ++ ) {
		site[column] = dataset->alphabet[getSiteSymbol(config, dataset, sequence,
							   offset, strand, column)];
	}
	return site;
}

// Write discovered motif sites in FASTA format
void writeMotifFASTA( const Config *config,
		      const Dataset *dataset,
		      const vector<PipelineResult> &pipelineResults,
		      const vector<OutputMotif> &outputMotifs ) {
	string filename = config->outputPrefix + ".fasta";
	ofstream outputFile(filename);
	if ( outputFile.is_open() == false ) {
		printf( "Unable to create output file: %s\n", filename.c_str() );
		exit(1);
	}

	if ( config->outputMotifNum == 1 && outputMotifs.empty() == false ) {
		const OutputMotif &motif = outputMotifs[0];
		const vector<uint32_t> &offsets = outputOffsets(pipelineResults, motif);
		const vector<uint8_t> &strands = outputStrands(pipelineResults, motif);
		for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx ++ ) {
			if ( motif.refined && motif.sitePresent[seqIdx] == 0 ) continue;
			outputFile << ">" << dataset->names[seqIdx] << " offset=" << offsets[seqIdx];
			if ( config->alphabetMode == ALPHABET_DNA ) {
				outputFile << " strand=" << (strands[seqIdx] == STRAND_REVERSE ? '-' : '+');
			}
			outputFile << "\n";
			outputFile << buildSiteString(config, dataset, dataset->sequences[seqIdx],
						  offsets[seqIdx], strands[seqIdx]) << "\n";
		}
	} else {
		for ( size_t motifIdx = 0; motifIdx < outputMotifs.size(); motifIdx ++ ) {
			const OutputMotif &motif = outputMotifs[motifIdx];
			int pipelineIdx = motif.pipelineIdx;
			const vector<uint32_t> &offsets = outputOffsets(pipelineResults, motif);
			const vector<uint8_t> &strands = outputStrands(pipelineResults, motif);
			for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx ++ ) {
				if ( motif.refined && motif.sitePresent[seqIdx] == 0 ) continue;
				outputFile << ">MorbiusPlus_" << motifIdx + 1 << "|" << seqIdx
					   << " motif_rank=" << motifIdx + 1
					   << " pipeline=" << pipelineIdx
					   << " offset=" << offsets[seqIdx];
				if ( config->alphabetMode == ALPHABET_DNA ) {
					outputFile << " strand=" << (strands[seqIdx] == STRAND_REVERSE ? '-' : '+');
				}
				outputFile << " source=" << dataset->names[seqIdx] << "\n";
				outputFile << buildSiteString(config, dataset, dataset->sequences[seqIdx],
							  offsets[seqIdx], strands[seqIdx]) << "\n";
			}
		}
	}
	closeOutput(outputFile, filename);
}

// Write sequence offsets
void writeOffsets( const Config *config,
		   const Dataset *dataset,
		   const vector<PipelineResult> &pipelineResults,
		   const vector<OutputMotif> &outputMotifs ) {
	string filename = config->outputPrefix + ".offsets.tsv";
	ofstream outputFile(filename);
	if ( outputFile.is_open() == false ) {
		printf( "Unable to create output file: %s\n", filename.c_str() );
		exit(1);
	}

	if ( config->outputMotifNum == 1 ) {
		outputFile << "SequenceIdx\tSequenceName\tOffset";
		if ( config->alphabetMode == ALPHABET_DNA ) outputFile << "\tStrand";
		outputFile << "\tMotif\n";
		if ( outputMotifs.empty() ) {
			closeOutput(outputFile, filename);
			return;
		}
		const OutputMotif &motif = outputMotifs[0];
		const vector<uint32_t> &offsets = outputOffsets(pipelineResults, motif);
		const vector<uint8_t> &strands = outputStrands(pipelineResults, motif);
		for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx ++ ) {
			if ( motif.refined && motif.sitePresent[seqIdx] == 0 ) continue;
			outputFile << seqIdx << "\t"
				   << dataset->names[seqIdx] << "\t"
				   << offsets[seqIdx] << "\t";
			if ( config->alphabetMode == ALPHABET_DNA ) {
				outputFile << (strands[seqIdx] == STRAND_REVERSE ? '-' : '+') << "\t";
			}
			outputFile << buildSiteString(config, dataset, dataset->sequences[seqIdx],
						  offsets[seqIdx], strands[seqIdx]) << "\n";
		}
	} else {
		outputFile << "MotifRank\tPipelineIdx\tSequenceIdx\tSequenceName\tOffset";
		if ( config->alphabetMode == ALPHABET_DNA ) outputFile << "\tStrand";
		outputFile << "\tMotif\n";
		for ( size_t motifIdx = 0; motifIdx < outputMotifs.size(); motifIdx ++ ) {
			const OutputMotif &motif = outputMotifs[motifIdx];
			int pipelineIdx = motif.pipelineIdx;
			const vector<uint32_t> &offsets = outputOffsets(pipelineResults, motif);
			const vector<uint8_t> &strands = outputStrands(pipelineResults, motif);
			for ( size_t seqIdx = 0; seqIdx < dataset->sequences.size(); seqIdx ++ ) {
				if ( motif.refined && motif.sitePresent[seqIdx] == 0 ) continue;
				outputFile << motifIdx + 1 << "\t"
					   << pipelineIdx << "\t"
					   << seqIdx << "\t"
					   << dataset->names[seqIdx] << "\t"
					   << offsets[seqIdx] << "\t";
				if ( config->alphabetMode == ALPHABET_DNA ) {
					outputFile << (strands[seqIdx] == STRAND_REVERSE ? '-' : '+') << "\t";
				}
				outputFile << buildSiteString(config, dataset, dataset->sequences[seqIdx],
							  offsets[seqIdx], strands[seqIdx]) << "\n";
			}
		}
	}
	closeOutput(outputFile, filename);
}

// Write a PWM table
void writePWM( const Config *config,
	       const Dataset *dataset,
	       const vector<OutputMotif> &outputMotifs ) {
	string filename = config->outputPrefix + ".pwm.tsv";
	ofstream outputFile(filename);
	if ( outputFile.is_open() == false ) {
		printf( "Unable to create output file: %s\n", filename.c_str() );
		exit(1);
	}
	if ( config->controlFilename.empty() == false ) outputFile << setprecision(17);

	if ( config->outputMotifNum == 1 ) {
		outputFile << "Position";
		for ( int symbol = 0; symbol < dataset->alphabetSize; symbol ++ ) {
			outputFile << "\t" << dataset->alphabet[symbol];
		}
		outputFile << "\n";

		if ( outputMotifs.empty() ) {
			closeOutput(outputFile, filename);
			return;
		}
		double denominator = (double)outputSiteNum(dataset, outputMotifs[0]) +
				     (double)dataset->alphabetSize * PSEUDOCOUNT;
		const vector<uint32_t> &count = outputMotifs[0].count;
		for ( size_t column = 0; column < config->motifLength; column ++ ) {
			outputFile << column;
			for ( int symbol = 0; symbol < dataset->alphabetSize; symbol ++ ) {
				double probability =
					((double)count[(size_t)symbol * config->motifLength + column] + PSEUDOCOUNT) /
					denominator;
				outputFile << "\t" << probability;
			}
			outputFile << "\n";
		}
	} else {
		outputFile << "MotifRank\tPipelineIdx\tPosition";
		for ( int symbol = 0; symbol < dataset->alphabetSize; symbol ++ ) {
			outputFile << "\t" << dataset->alphabet[symbol];
		}
		outputFile << "\n";

		for ( size_t motifIdx = 0; motifIdx < outputMotifs.size(); motifIdx ++ ) {
			double denominator = (double)outputSiteNum(dataset, outputMotifs[motifIdx]) +
					     (double)dataset->alphabetSize * PSEUDOCOUNT;
			const vector<uint32_t> &count = outputMotifs[motifIdx].count;
			for ( size_t column = 0; column < config->motifLength; column ++ ) {
				outputFile << motifIdx + 1 << "\t"
					   << outputMotifs[motifIdx].pipelineIdx << "\t"
					   << column;
				for ( int symbol = 0; symbol < dataset->alphabetSize; symbol ++ ) {
					double probability =
						((double)count[(size_t)symbol * config->motifLength + column] + PSEUDOCOUNT) /
						denominator;
					outputFile << "\t" << probability;
				}
				outputFile << "\n";
			}
		}
	}
	closeOutput(outputFile, filename);
}

// Write a MEME-format motif
void writeMEME( const Config *config,
		const Dataset *dataset,
		const vector<OutputMotif> &outputMotifs ) {
	string filename = config->outputPrefix + ".meme";
	ofstream outputFile(filename);
	if ( outputFile.is_open() == false ) {
		printf( "Unable to create output file: %s\n", filename.c_str() );
		exit(1);
	}
	if ( config->controlFilename.empty() == false ) outputFile << setprecision(17);

	outputFile << "MEME version 4\n\n";
	outputFile << "ALPHABET= " << dataset->alphabet << "\n\n";
	if ( config->alphabetMode == ALPHABET_DNA ) outputFile << "strands: +\n\n";
	outputFile << "Background letter frequencies\n";
	for ( int symbol = 0; symbol < dataset->alphabetSize; symbol ++ ) {
		outputFile << dataset->alphabet[symbol] << " " << 1.0 / (double)dataset->alphabetSize;
		if ( symbol + 1 < dataset->alphabetSize ) outputFile << " ";
	}
	outputFile << "\n\n";
	for ( size_t motifIdx = 0; motifIdx < outputMotifs.size(); motifIdx ++ ) {
		size_t siteNum = outputSiteNum(dataset, outputMotifs[motifIdx]);
		double denominator = (double)siteNum + (double)dataset->alphabetSize * PSEUDOCOUNT;
		if ( config->outputMotifNum == 1 ) outputFile << "MOTIF MorbiusPlus\n";
		else outputFile << "MOTIF MorbiusPlus_" << motifIdx + 1 << "\n";
		outputFile << "letter-probability matrix: alength= " << dataset->alphabetSize
			   << " w= " << config->motifLength
			   << " nsites= " << siteNum
			   << " E= 0\n";

		const vector<uint32_t> &count = outputMotifs[motifIdx].count;
		for ( size_t column = 0; column < config->motifLength; column ++ ) {
			for ( int symbol = 0; symbol < dataset->alphabetSize; symbol ++ ) {
				double probability =
					((double)count[(size_t)symbol * config->motifLength + column] + PSEUDOCOUNT) /
					denominator;
				outputFile << probability;
				if ( symbol + 1 < dataset->alphabetSize ) outputFile << " ";
			}
			outputFile << "\n";
		}
		if ( motifIdx + 1 < outputMotifs.size() ) outputFile << "\n";
	}
	closeOutput(outputFile, filename);
}

// Write a result summary
void writeSummary( const Config *config,
		   const Dataset *dataset,
		   const vector<SeedModel> &seedModels,
		   const vector<PipelineResult> &pipelineResults,
		   const vector<OutputMotif> &outputMotifs,
		   double elapsedTime ) {
	string filename = config->outputPrefix + ".summary.txt";
	ofstream outputFile(filename);
	if ( outputFile.is_open() == false ) {
		printf( "Unable to create output file: %s\n", filename.c_str() );
		exit(1);
	}

	outputFile << "Input File              : " << config->inputFilename << "\n";
	outputFile << "Alphabet                : " << (config->alphabetMode == ALPHABET_DNA ? "DNA" : "Protein") << "\n";
	outputFile << "Sequence Number          : " << dataset->sequences.size() << "\n";
	outputFile << "Sequence Length          : " << dataset->sequenceLength << "\n";
	outputFile << "Motif Length             : " << config->motifLength << "\n";
	if ( outputMotifs.empty() == false ) {
		const SeedModel *seedModel = &seedModels[outputMotifs[0].pipelineIdx];
		outputFile << "SampleNum                : " << seedModel->sampleNum << "\n";
		outputFile << "Seed Pipeline            : " << outputMotifs[0].pipelineIdx << "\n";
		outputFile << "Seed Valid               : " << (seedModel->valid ? "Yes" : "No") << "\n";
		if ( seedModel->valid ) {
			outputFile << "Selected Seed            : " << seedModel->seedString << "\n";
			outputFile << "Seed Support             : " << seedModel->seedSupport << "\n";
			outputFile << "Expected Seed Support    : " << seedModel->seedExpectedSupport << "\n";
			outputFile << "Seed Rank                : " << seedModel->seedRank << "\n";
			outputFile << "Anchor Offset            : " << seedModel->anchorOffset << "\n";
			outputFile << "Anchor Rank              : " << seedModel->anchorRank << "\n";
		}
	}
	outputFile << "Pipeline Number          : " << NUMPIPELINE << "\n";
	if ( config->controlFilename.empty() == false ) {
		outputFile << "Refinement               : Separate post-Gibbs Primary/Control enrichment\n";
		outputFile << "Best Score Meaning       : Original Gibbs agreement score (all Primary sequences)\n";
		outputFile << "Elapsed Time Includes    : Input, seed initialization, Gibbs, refinement, and result output before summary\n";
	}
	if ( config->outputMotifNum == 1 && outputMotifs.empty() == false ) {
		const OutputMotif &outputMotif = outputMotifs[0];
		const PipelineResult &bestResult = pipelineResults[outputMotif.pipelineIdx];
		outputFile << "Selected Pipeline        : " << outputMotif.pipelineIdx << "\n";
		outputFile << "Best Score               : " << bestResult.bestScore << "\n";
		outputFile << "Normalized Best Score    : "
			   << calculateNormalizedScore(config, dataset, bestResult.bestScore) << "\n";
		outputFile << "Score Threshold          : " << config->scoreThreshold << "\n";
		outputFile << "Pipeline Update Number   : " << bestResult.updateNum << "\n";
		outputFile << "Threshold Reached        : " << (bestResult.thresholdReached ? "Yes" : "No") << "\n";
		outputFile << "Consensus Subsequence    : " << outputMotif.consensus << "\n";
		if ( outputMotif.refined ) outputFile << "Refined Site Number      : " << outputMotif.siteNum << "\n";
	} else {
		outputFile << "Score Threshold          : " << config->scoreThreshold << "\n";
		outputFile << "Requested Motif Number   : " << config->outputMotifNum << "\n";
		outputFile << "Reported Motif Number    : " << outputMotifs.size() << "\n";
		for ( size_t motifIdx = 0; motifIdx < outputMotifs.size(); motifIdx ++ ) {
			const OutputMotif &outputMotif = outputMotifs[motifIdx];
			const PipelineResult &pipelineResult = pipelineResults[outputMotif.pipelineIdx];
			outputFile << "\n[Motif " << motifIdx + 1 << "]\n";
			outputFile << "Motif ID                 : MorbiusPlus_" << motifIdx + 1 << "\n";
			outputFile << "Selected Pipeline        : " << outputMotif.pipelineIdx << "\n";
			outputFile << "Best Score               : " << pipelineResult.bestScore << "\n";
			outputFile << "Normalized Best Score    : "
				   << calculateNormalizedScore(config, dataset, pipelineResult.bestScore) << "\n";
			outputFile << "Pipeline Update Number   : " << pipelineResult.updateNum << "\n";
			outputFile << "Threshold Reached        : "
				   << (pipelineResult.thresholdReached ? "Yes" : "No") << "\n";
			outputFile << "Consensus Subsequence    : " << outputMotif.consensus << "\n";
			if ( outputMotif.refined ) outputFile << "Refined Site Number      : " << outputMotif.siteNum << "\n";
		}
	}
	outputFile << "Elapsed Time             : " << elapsedTime << " seconds\n";
	if ( config->alphabetMode == ALPHABET_DNA ) {
		outputFile << "\nStrand Search            : Forward only\n";
		outputFile << "Initial Strand           : + (all sequences and pipelines)\n";
		outputFile << "Offset Coordinates       : Original sequence, zero-based leftmost position\n";
		outputFile << "Maximum Pipeline Updates : " << config->maxUpdateNum << "\n";
		outputFile << "\n[All Pipelines]\n";
		outputFile << "Pipeline\tBestScore\tUpdates\tThresholdReached\tTerminationReason\n";
		for ( size_t pipelineIdx = 0; pipelineIdx < pipelineResults.size(); pipelineIdx ++ ) {
			const PipelineResult &pipelineResult = pipelineResults[pipelineIdx];
			outputFile << pipelineIdx << "\t" << pipelineResult.bestScore << "\t"
				   << pipelineResult.updateNum << "\t"
				   << (pipelineResult.thresholdReached ? "Yes" : "No") << "\t"
				   << (pipelineResult.thresholdReached ? "threshold" : "max_updates") << "\n";
		}
	}
	closeOutput(outputFile, filename);
}

// Print the final result
void printResult( const Config *config,
		  const Dataset *dataset,
		  const vector<SeedModel> &seedModels,
		  const vector<PipelineResult> &pipelineResults,
		  const vector<OutputMotif> &outputMotifs,
		  double elapsedTime ) {
	printf( "---------------------------------------------------------------------\n" );
	printf( "MORBIUS+ RESULT\n" );
	printf( "---------------------------------------------------------------------\n" );
	printf( "The Number of Sequence : %lu\n", (unsigned long)dataset->sequences.size() );
	printf( "The Length of Sequence : %lu\n", (unsigned long)dataset->sequenceLength );
	printf( "The Length of Motif    : %lu\n", (unsigned long)config->motifLength );
	if ( outputMotifs.empty() == false ) {
		const SeedModel *seedModel = &seedModels[outputMotifs[0].pipelineIdx];
		printf( "SampleNum              : %lu\n", (unsigned long)seedModel->sampleNum );
		printf( "Seed Pipeline          : %d\n", outputMotifs[0].pipelineIdx );
		if ( seedModel->valid ) {
			printf( "Selected Seed          : %s\n", seedModel->seedString.c_str() );
			printf( "Anchor Offset          : %lu\n", (unsigned long)seedModel->anchorOffset );
		} else {
			printf( "Selected Seed          : None\n" );
		}
	}
	if ( config->controlFilename.empty() == false ) {
		printf( "Refinement             : Separate post-Gibbs Primary/Control enrichment\n" );
		printf( "Best Score Meaning     : Original Gibbs agreement score (all Primary sequences)\n" );
	}
	if ( config->outputMotifNum == 1 && outputMotifs.empty() == false ) {
		const OutputMotif &outputMotif = outputMotifs[0];
		const PipelineResult &bestResult = pipelineResults[outputMotif.pipelineIdx];
		printf( "Selected Pipeline      : %d\n", outputMotif.pipelineIdx );
		printf( "Best Score             : %lu\n", (unsigned long)bestResult.bestScore );
		printf( "Normalized Best Score  : %.8f\n",
			calculateNormalizedScore(config, dataset, bestResult.bestScore) );
		printf( "Pipeline Updates       : %lu\n", (unsigned long)bestResult.updateNum );
		printf( "Consensus Subsequence  : %s\n", outputMotif.consensus.c_str() );
		if ( outputMotif.refined ) printf( "Refined Site Number    : %lu\n", (unsigned long)outputMotif.siteNum );
	} else {
		printf( "Requested Motifs       : %lu\n", (unsigned long)config->outputMotifNum );
		printf( "Reported Motifs        : %lu\n", (unsigned long)outputMotifs.size() );
		for ( size_t motifIdx = 0; motifIdx < outputMotifs.size(); motifIdx ++ ) {
			const OutputMotif &outputMotif = outputMotifs[motifIdx];
			const PipelineResult &pipelineResult = pipelineResults[outputMotif.pipelineIdx];
			printf( "[Motif %lu]\n", (unsigned long)motifIdx + 1 );
			printf( "Selected Pipeline      : %d\n", outputMotif.pipelineIdx );
			printf( "Best Score             : %lu\n", (unsigned long)pipelineResult.bestScore );
			printf( "Normalized Best Score  : %.8f\n",
				calculateNormalizedScore(config, dataset, pipelineResult.bestScore) );
			printf( "Pipeline Updates       : %lu\n", (unsigned long)pipelineResult.updateNum );
			printf( "Consensus Subsequence  : %s\n", outputMotif.consensus.c_str() );
			if ( outputMotif.refined ) printf( "Refined Site Number    : %lu\n", (unsigned long)outputMotif.siteNum );
		}
	}
	printf( "Elapsed Time           : %.8f\n", elapsedTime );
	printf( "---------------------------------------------------------------------\n" );
	fflush( stdout );
}
