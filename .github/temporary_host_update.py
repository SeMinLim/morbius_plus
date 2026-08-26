from pathlib import Path


def replace_once(text, old, new, label):
	count = text.count(old)
	if count != 1:
		raise RuntimeError(f"{label}: expected one match, found {count}")
	return text.replace(old, new, 1)


def replace_between(text, start, end, replacement, label):
	begin = text.find(start)
	if begin < 0:
		raise RuntimeError(f"{label}: start marker not found")
	finish = text.find(end, begin)
	if finish < 0:
		raise RuntimeError(f"{label}: end marker not found")
	return text[:begin] + replacement + text[finish:]


PARSE_RESULT = r'''// Parse one set of per-sequence result beats
void parseResultBeat( const uint8_t *beat,
			      vector<AccelWireResult> &results ) {
	results.assign(ACCELMAXPIPELINE, AccelWireResult());
	for ( int resultBeatIdx = 0; resultBeatIdx < ACCELRESULTBEATNUM; resultBeatIdx ++ ) {
		const uint8_t *resultBeat = beat + (size_t)resultBeatIdx * ACCELBEATBYTES;
		if ( readUInt32(resultBeat + 0) != ACCELRESULTMAGIC ||
		     readUInt16(resultBeat + 4) != ACCELPROTOCOLVERSION ||
		     resultBeat[6] != (uint8_t)resultBeatIdx ) {
			printf( "Invalid Morbius+ accelerator result header at result beat %d.\n",
				resultBeatIdx
			);
			exit(1);
		}
		for ( int resultSlot = 0;
		      resultSlot < ACCELPIPELINESPERRESULTBEAT;
		      resultSlot ++ ) {
			int pipelineIdx = resultBeatIdx * ACCELPIPELINESPERRESULTBEAT + resultSlot;
			size_t base = 64 + (size_t)resultSlot * 96;
			results[(size_t)pipelineIdx].newOffset = getBits(resultBeat, base + 0, 11);
			results[(size_t)pipelineIdx].bestUpdate = getBits(resultBeat, base + 11, 1) != 0;
			results[(size_t)pipelineIdx].terminated = getBits(resultBeat, base + 12, 1) != 0;
			results[(size_t)pipelineIdx].active = getBits(resultBeat, base + 13, 1) != 0;
			results[(size_t)pipelineIdx].bestScore = getBits(resultBeat, base + 16, 32);
			results[(size_t)pipelineIdx].updateNum = getBits(resultBeat, base + 48, 32);
		}
	}
}

'''

PARSE_SUMMARY = r'''// Parse one fixed-position summary beat
void parseSummaryBeat( const uint8_t *beat, AccelSummary *summary ) {
	if ( readUInt32(beat + 0) != ACCELRESULTMAGIC ||
	     readUInt16(beat + 4) != ACCELPROTOCOLVERSION ||
	     beat[6] != 0xffU ) {
		printf( "Invalid Morbius+ accelerator summary header.\n" );
		exit(1);
	}
	summary->processedNum = readUInt32(beat + 8);
	summary->allDone = (beat[12] & 0x1U) != 0;
	summary->bestPipelineIdx = beat[13];
	summary->cycleNum = readUInt32(beat + 16);
}
'''

SELECT_LOCAL_MODEL = r'''static uint32_t selectLocalCandidateModel( const uint32_t weight[ACCELSEGMENTSIZE],
					   uint32_t totalMass,
					   uint32_t randomFraction,
					   uint32_t validNum ) {
	uint64_t product = (uint64_t)totalMass * (uint64_t)randomFraction;
	uint32_t threshold = (uint32_t)(product >> 24);
	uint32_t cumulative = 0;
	for ( uint32_t localIdx = 0; localIdx < validNum; localIdx ++ ) {
		uint32_t nextCumulative = cumulative + weight[localIdx];
		if ( nextCumulative > threshold || localIdx + 1 == validNum ) return localIdx;
		cumulative = nextCumulative;
	}
	return 0;
}

'''

PACK_RESULT_MODEL = r'''static void packResultModel( uint8_t *output,
			     const AccelWireResult result[ACCELMAXPIPELINE],
			     bool allDone,
			     int bestPipelineIdx ) {
	for ( int resultBeatIdx = 0; resultBeatIdx < ACCELRESULTBEATNUM; resultBeatIdx ++ ) {
		uint8_t *word = output + (size_t)resultBeatIdx * ACCELBEATBYTES;
		writeUInt32Model(word + 0, ACCELRESULTMAGIC);
		writeUInt16Model(word + 4, ACCELPROTOCOLVERSION);
		word[6] = (uint8_t)resultBeatIdx;
		for ( int resultSlot = 0;
		      resultSlot < ACCELPIPELINESPERRESULTBEAT;
		      resultSlot ++ ) {
			int pipelineIdx = resultBeatIdx * ACCELPIPELINESPERRESULTBEAT + resultSlot;
			size_t base = 64 + (size_t)resultSlot * 96;
			setBitsModel(word, base + 0, result[pipelineIdx].newOffset, 11);
			setBitsModel(word, base + 11, result[pipelineIdx].bestUpdate ? 1 : 0, 1);
			setBitsModel(word, base + 12, result[pipelineIdx].terminated ? 1 : 0, 1);
			setBitsModel(word, base + 13, result[pipelineIdx].active ? 1 : 0, 1);
			setBitsModel(word, base + 16, result[pipelineIdx].bestScore, 32);
			setBitsModel(word, base + 48, result[pipelineIdx].updateNum, 32);
		}
	}
	(void)allDone;
	(void)bestPipelineIdx;
}

'''

PACK_SUMMARY_MODEL = r'''static void packSummaryModel( uint8_t *word,
			      uint32_t processedNum,
			      bool allDone,
			      int bestPipelineIdx ) {
	writeUInt32Model(word + 0, ACCELRESULTMAGIC);
	writeUInt16Model(word + 4, ACCELPROTOCOLVERSION);
	word[6] = 0xffU;
	writeUInt32Model(word + 8, processedNum);
	word[12] = allDone ? 1 : 0;
	word[13] = (uint8_t)bestPipelineIdx;
	writeUInt32Model(word + 16, 0);
}

'''


def update_header(path):
	text = path.read_text()
	text = replace_once(text, "#define DEFAULTNUMPIPELINE 4", "#define DEFAULTNUMPIPELINE 16", "default pipelines")
	text = replace_once(text, "#define ACCELMAXPIPELINE 4", "#define ACCELMAXPIPELINE 16", "physical pipelines")
	text = replace_once(text, "#define ACCELSEGMENTSIZE 8", "#define ACCELSEGMENTSIZE 128", "profiler segment size")
	text = replace_once(text, "#define ACCELPROTOCOLVERSION 1", "#define ACCELPROTOCOLVERSION 2", "protocol version")
	text = replace_once(
		text,
		"#define ACCELBEATBYTES 64",
		"#define ACCELBEATBYTES 64\n"
		"#define ACCELPIPELINESPERRESULTBEAT 4\n"
		"#define ACCELRESULTBEATNUM (ACCELMAXPIPELINE / ACCELPIPELINESPERRESULTBEAT)\n"
		"#define ACCELSUMMARYBEATNUM 1",
		"result protocol constants",
	)
	path.write_text(text)


def update_protocol(path):
	text = path.read_text()
	text = replace_between(text, "// Parse one per-sequence result beat", "// Parse one fixed-position summary beat", PARSE_RESULT, "result parser")
	begin = text.find("// Parse one fixed-position summary beat")
	if begin < 0:
		raise RuntimeError("summary parser marker not found")
	text = text[:begin] + PARSE_SUMMARY
	path.write_text(text)


def update_orchestrator(path):
	text = path.read_text()
	text = replace_once(
		text,
		"if ( executor(input, 2 * ACCELBEATBYTES, output, executorContext) != 0 ) {",
		"if ( executor(input,\n"
		"\t\t     (ACCELRESULTBEATNUM + ACCELSUMMARYBEATNUM) * ACCELBEATBYTES,\n"
		"\t\t     output,\n"
		"\t\t     executorContext) != 0 ) {",
		"bootstrap output size",
	)
	text = replace_once(
		text,
		"parseSummaryBeat(output.data() + ACCELBEATBYTES, &summary);",
		"parseSummaryBeat(output.data() + ACCELRESULTBEATNUM * ACCELBEATBYTES, &summary);",
		"bootstrap summary offset",
	)
	text = replace_once(
		text,
		"size_t outputSize = (batchNum + 1) * ACCELBEATBYTES;",
		"size_t outputSize = (batchNum * ACCELRESULTBEATNUM + ACCELSUMMARYBEATNUM) *\n"
		"\t\t\t    ACCELBEATBYTES;",
		"batch output size",
	)
	text = replace_once(
		text,
		"parseSummaryBeat(output.data() + batchNum * ACCELBEATBYTES, &summary);",
		"parseSummaryBeat(output.data() + batchNum * ACCELRESULTBEATNUM * ACCELBEATBYTES,\n"
		"\t\t\t &summary);",
		"batch summary offset",
	)
	text = replace_once(
		text,
		"parseResultBeat(output.data() + (size_t)itemIdx * ACCELBEATBYTES, wireResults);",
		"parseResultBeat(output.data() +\n"
		"\t\t\t\t(size_t)itemIdx * ACCELRESULTBEATNUM * ACCELBEATBYTES,\n"
		"\t\t\t\twireResults);",
		"batch result offset",
	)
	path.write_text(text)


def update_model(path):
	text = path.read_text()
	text = replace_between(text, "static uint32_t selectLocalCandidateModel", "static uint32_t sampleCandidateModel", SELECT_LOCAL_MODEL, "model candidate selector")
	text = text.replace("uint32_t logProb[ACCELSEGMENTSIZE] = {0, 0, 0, 0, 0, 0, 0, 0};", "uint32_t logProb[ACCELSEGMENTSIZE] = {0};")
	text = text.replace("uint32_t weight[ACCELSEGMENTSIZE] = {0, 0, 0, 0, 0, 0, 0, 0};", "uint32_t weight[ACCELSEGMENTSIZE] = {0};")
	text = replace_between(text, "static void packResultModel", "static void packSummaryModel", PACK_RESULT_MODEL, "model result packer")
	text = replace_between(text, "static void packSummaryModel", "void resetModelKernel", PACK_SUMMARY_MODEL, "model summary packer")
	text = replace_once(
		text,
		"if ( outputSize < ((size_t)batchSize + 1) * ACCELBEATBYTES ) {",
		"if ( outputSize < ((size_t)batchSize * ACCELRESULTBEATNUM +\n"
		"\t\t\t   ACCELSUMMARYBEATNUM) * ACCELBEATBYTES ) {",
		"model output size validation",
	)
	old_call = '''\t\t\tpackResultModel(output.data() + (size_t)itemIdx * ACCELBEATBYTES,
\t\t\t\t\t(uint16_t)itemIdx,
\t\t\t\t\tresult,
\t\t\t\t\tallDone,
\t\t\t\t\tbestPipelineIdx);'''
	new_call = '''\t\t\tpackResultModel(output.data() +
\t\t\t\t\t(size_t)itemIdx * ACCELRESULTBEATNUM * ACCELBEATBYTES,
\t\t\t\t\tresult,
\t\t\t\t\tallDone,
\t\t\t\t\tbestPipelineIdx);'''
	text = replace_once(text, old_call, new_call, "model batch result call")
	text = replace_once(
		text,
		"packResultModel(output.data(), 0, result, allDone, bestPipelineIdx);",
		"packResultModel(output.data(), result, allDone, bestPipelineIdx);",
		"model bootstrap result call",
	)
	text = replace_once(
		text,
		"packSummaryModel(output.data() + (size_t)batchSize * ACCELBEATBYTES,",
		"packSummaryModel(output.data() +\n"
		"\t\t\t (size_t)batchSize * ACCELRESULTBEATNUM * ACCELBEATBYTES,",
		"model summary offset",
	)
	path.write_text(text)


def update_protocol_test(path):
	text = path.read_text()
	text = replace_once(text, "\tconfig->numPipeline = 1;", "\tconfig->numPipeline = ACCELMAXPIPELINE;", "test pipeline count")
	text = replace_once(text, "\t\tstate.terminated = pipelineIdx != 0;", "\t\tstate.terminated = false;", "test active states")
	text = text.replace("executeModelKernel(input, 2 * ACCELBEATBYTES, output, NULL)", "executeModelKernel(input,\n\t\t\t       (ACCELRESULTBEATNUM + ACCELSUMMARYBEATNUM) * ACCELBEATBYTES,\n\t\t\t       output,\n\t\t\t       NULL)")
	old_boot = '''\tif ( results[0].newOffset != 12 || results[0].bestScore != 36 ||
\t     results[0].updateNum != 1 || results[0].bestUpdate == false ||
\t     results[0].terminated ) {'''
	new_boot = '''\tif ( results.size() != ACCELMAXPIPELINE || results[0].newOffset > 12 ||
\t     results[0].bestScore < 32 || results[0].updateNum != 1 ||
\t     results[0].terminated ) {'''
	text = replace_once(text, old_boot, new_boot, "protocol bootstrap assertion")
	old_update = '''\tif ( results[0].newOffset != 4 || results[0].bestScore != 36 ||
\t     results[0].updateNum != 2 || results[0].bestUpdate ||
\t     results[0].terminated == false ) {'''
	new_update = '''\tif ( results.size() != ACCELMAXPIPELINE || results[0].newOffset > 12 ||
\t     results[0].bestScore < 32 || results[0].updateNum != 2 ||
\t     results[0].terminated == false ) {'''
	text = replace_once(text, old_update, new_update, "protocol update assertion")
	path.write_text(text)


def update_host_readme(path):
	text = path.read_text().replace("--pipelines 4", "--pipelines 16")
	path.write_text(text)


def update_host_tree(root):
	update_header(root / "MorbiusPlus.h")
	update_protocol(root / "AcceleratorProtocol.cpp")
	update_orchestrator(root / "HostOrchestrator.cpp")
	update_model(root / "AcceleratorModel.cpp")
	update_protocol_test(root / "test/test_protocol.cpp")
	update_host_readme(root / "README.md")
