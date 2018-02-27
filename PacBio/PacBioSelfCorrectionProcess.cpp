///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioSelfCorrectionProcess.cpp - Self-correction using FM-index walk for PacBio reads
//
#include "PacBioSelfCorrectionProcess.h"
#include "SAIPBSelfCTree.h"
#include "CorrectionThresholds.h"
#include "HashMap.h"
#include "SAIPBHybridCTree.h"
#include "LongReadOverlap.h"
#include "Timer.h"
#include "KmerThresholdTable.h"
#include "Util.h"
#include "Kmer.h"
#include "Alphabet.h"
#include <iomanip>
#include <time.h>
#include <algorithm>
#include <memory>

#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)

// PacBio Self Correction by Ya and YTH, v20151202.
// 1. Identify highly-accurate seeds within PacBio reads
// 2. For each pair of seeds, perform kmer extension using local kmer frequency collected by FM-index extension
PacBioSelfCorrectionResult PacBioSelfCorrectionProcess::process(const SequenceWorkItem& workItem)
{
	PacBioSelfCorrectionResult result;
    result.readid = workItem.read.id;
	std::string readSeq = workItem.read.seq.toString();
	const size_t readSeqLen = readSeq.length();
	SeedVector seedVec, pieceVec;

	for(auto& iter : m_params.kmerSet)
		Kmer::kmerMap[iter] = std::unique_ptr<Kmer[]>(new Kmer[readSeqLen]);
	searchSeedsWithHybridKmers(readSeq, seedVec, result);
    initCorrect(readSeq, seedVec, pieceVec, result);
	Kmer::kmerMap.clear();

	result.merge = !pieceVec.empty();
	result.totalReadsLen = readSeq.length();
	for(const auto& iter : pieceVec)
		result.correctedStrs.push_back(iter.seedStr);
	return result;
}

// Search seeds with static and dynamic kmers. Noted by KuanWeiLee 20171027
void PacBioSelfCorrectionProcess::searchSeedsWithHybridKmers(const std::string& readSeq, SeedVector& seedVec, PacBioSelfCorrectionResult& result)
{
	const size_t readSeqLen = readSeq.length();
	int staticKmerSize = m_params.startKmerLength;
	if((int)readSeqLen < staticKmerSize) return;

    Timer* seedTimer = new Timer("Seed Time", true);
	int attribute[readSeqLen];
	getSeqAttribute(readSeq, attribute);
	if(m_params.Manual) std::fill_n(attribute, readSeqLen, m_params.mode);

	//Search seeds; slide through the read sequence with hybrid-kmers. Noted by KuanWeiLee
	//kmer[Start/Move]Pos indicate the starting/moving position of the static-kmer.
	float** table = KmerThresholdTable::m_table;

	for(size_t kmerStartPos = 0; kmerStartPos < readSeqLen; kmerStartPos++)
	{
		int kmerStartType = attribute[kmerStartPos];
		staticKmerSize -= m_params.kmerOffset[kmerStartType];
		bool isSeed = false, isRepeat = false;
		Kmer dynamicKmer = Kmer::kmerMap[staticKmerSize][kmerStartPos];
	//	Kmer dynamicKmer(&m_params.indices, readSeq, kmerStartPos, staticKmerSize);
		int maxFixedMerFreq = dynamicKmer.getFreq();
		size_t seedStartPos = kmerStartPos;
		for(size_t kmerMovePos = kmerStartPos; kmerMovePos < readSeqLen; kmerMovePos++)
		{
			int kmerMoveType = attribute[kmerMovePos];
			const Kmer& staticKmer = Kmer::kmerMap[staticKmerSize][kmerMovePos];
		//	const Kmer staticKmer(&m_params.indices, readSeq, kmerMovePos, staticKmerSize);
			if(isSeed)
			{
				char b = readSeq[(kmerMovePos + staticKmerSize - 1)];
				dynamicKmer.expand(b);
			}
			float dynamicThreshold = table[kmerStartType][dynamicKmer.getSize()];
			float staticThreshold = table[kmerMoveType][staticKmerSize];
			float repeatThreshold = staticThreshold * (5 - ((kmerMoveType >> 1) << 2));
			bool isOnRepeat = (staticKmer.getFreq() >= repeatThreshold);
			float freqDiff = (float)staticKmer.getFreq()/maxFixedMerFreq;
			//Gerneral seed extension strategy.
			if	(
				   dynamicKmer.getSize() > m_params.kmerLengthUpperBound		//1.over length
				|| !dynamicKmer.isValid()										//2.dynamic frequency(1)
				|| dynamicKmer.getFreq() < dynamicThreshold						//2.dynamic frequency(2)
				|| staticKmer.getFreq() < staticThreshold						//3.static frequency
				)
			{
				if(isSeed && !staticKmer.getPseudo()) dynamicKmer.shrink(1);
				break;
			}
			//Kmer Hitchhike strategy.
			int isGiantRepeat = ((kmerStartType >> 1) & (kmerMoveType >> 1)) + 1;
			if(isRepeat && freqDiff < (m_params.hhRatio/isGiantRepeat))			//4.hitchhiking kmer(1) (HIGH-->LOW)
			{
				dynamicKmer.shrink(1);
				kmerStartPos++;
				break;
			}
			else if(isOnRepeat && freqDiff > (isGiantRepeat/m_params.hhRatio))	//4.hitchhiking kmer(2) (LOW-->HIGH)
			{
				isSeed = false;
				kmerStartPos = kmerMovePos - 1;
				break;
			}
			isSeed = true;
			kmerStartPos = seedStartPos + dynamicKmer.getSize() - 1;
			isRepeat = isRepeat || isOnRepeat;
			maxFixedMerFreq = MAX(maxFixedMerFreq, staticKmer.getFreq());
		}
		//Low Complexity strategy.
		if(isSeed && !dynamicKmer.isLowComplexity())
		{
			SeedFeature newSeed(dynamicKmer.getWord(), seedStartPos, maxFixedMerFreq, isRepeat, staticKmerSize, m_params.PBcoverage);
			newSeed.estimateBestKmerSize(m_params.indices);
			seedVec.push_back(newSeed);
		}
		staticKmerSize += m_params.kmerOffset[kmerStartType];
	}

	//Seed Hitchhike strategy.
	seedVec = removeHitchhikingSeeds(seedVec, attribute, result);
	if(m_params.DebugSeed)
	{
		std::ostream* pSeedWriter = createWriter(m_params.directory + "seed/" + result.readid + ".seed");
		write(*pSeedWriter, seedVec);
		delete pSeedWriter;
	}

	result.totalSeedNum = seedVec.size();
	result.Timer_Seed = seedTimer->getElapsedWallTime();
    delete seedTimer;
}
//Sequence attribute is set dynamically using a sliding fixed-mer on each position of the sequence.
//Noted by KuanWeiLee 20180118
void PacBioSelfCorrectionProcess::getSeqAttribute(const std::string& seq, int* const attribute)
{
	const size_t seqlen = seq.length();
	std::fill_n(attribute, seqlen, 1);

	int range = 300;
	int x = m_params.PBcoverage;
	int y = m_params.scanKmerLength;
	const int ksize = m_params.scanKmerLength;
//	float lowcov = KmerThresholdTable::calculate(0, x, y);
//	float unique = KmerThresholdTable::calculate(1, x, y);
	float repeat = KmerThresholdTable::calculate(2, x ,y);
	int front = 0, fear = -1;
	int leftmost = (seqlen - 1), rightmost = 0, repeatcount = 0;
	std::map<int, int> set;

	for(size_t pos = 0; pos < seqlen; pos++)
	{
		int left = pos - (range >> 1);
		int right = pos + (range >> 1);
		left = MAX(left, 0);
		right = MIN(right, (int)(seqlen - 1));
		while(fear < right)
		{
			fear++;
			Kmer* prev = nullptr;
			for(auto& iter : m_params.kmerSet)
			{
				Kmer::kmerMap[iter][fear] = Kmer(&m_params.indices, seq, fear, iter, prev);
				prev = Kmer::kmerMap[iter].get() + fear;
			}
			const Kmer& inKmer = Kmer::kmerMap[ksize][fear];
		//	Kmer inKmer(&m_params.indices, seq, ++fear, ksize);
			int freq = inKmer.isLowComplexity() ? 0 : inKmer.getFreq();
		//	int freq = inKmer.getFreq();
			int type;
			if(freq == 0) type = -1;
			else if(freq >= repeat) type = 2;
			else type = 1;
			set[type]++;
		}
		while(front < left)
		{
			const Kmer& outKmer = Kmer::kmerMap[ksize][front];
			front++;
		//	Kmer outKmer(&m_params.indices, seq, front++, ksize);
			int freq = outKmer.isLowComplexity() ? 0 : outKmer.getFreq();
		//	int freq = outKmer.getFreq();
			int type;
			if(freq == 0) type = -1;
			else if(freq >= repeat) type = 2;
			else type = 1;
			set[type]--;
		}
		int size = (right - left + 1) - set[-1];
		float ratio = (float)set[2]/size + 0.0005;
		if(ratio >= 0.02)
		{
			attribute[pos] = 2;
			repeatcount++;
			leftmost = MIN(leftmost, (int)pos);
			rightmost = MAX(rightmost, (int)pos);
		}
	}

	if((float)repeatcount/seqlen >= 0.5 && (float)(leftmost + (seqlen -rightmost))/seqlen <= 0.1)
		std::fill_n(attribute, seqlen, 2);
}

//Kmer & Seed Hitchhike strategy would maitain seed-correctness,
//once the sequence is stuck between the ambiguity from uniqu to repeat mode.
//Noted by KuanWeiLee 20180106
PacBioSelfCorrectionProcess::SeedVector PacBioSelfCorrectionProcess::removeHitchhikingSeeds(SeedVector initSeedVec, int const *attribute, PacBioSelfCorrectionResult& result)
{
	if(initSeedVec.size() < 2) return initSeedVec;
	int x = m_params.PBcoverage;
	int y = m_params.startKmerLength;
	float overfrequency = KmerThresholdTable::calculate(2, x ,y) * 5;
	for(SeedVector::iterator iterQuery = initSeedVec.begin(); (iterQuery + 1) != initSeedVec.end(); iterQuery++)
	{
		SeedFeature& query = *iterQuery;
		SeedVector::iterator iterTarget = iterQuery + 1;
		//if(query.isHitchhiked) continue;
		int queryType = attribute[query.seedStartPos];
		if(queryType == 2 && query.maxFixedMerFreq >= overfrequency) continue;
		for(; iterTarget != initSeedVec.end(); iterTarget++)
		{
			SeedFeature& target = *iterTarget;
			//if(target.isHitchhiked) continue;
			int	targetType = attribute[target.seedStartPos];
		//	int isBoundary = (queryType >> 1) ^ (targetType >> 1);
		//	if((int)(target.seedStartPos - query.seedEndPos) > (m_params.repeatDistance >> isBoundary)) break;
			if((int)(target.seedStartPos - query.seedEndPos) > m_params.repeatDistance) break;
			if(targetType == 2 && target.maxFixedMerFreq >= overfrequency) continue;
			float freqDiff = (float)target.maxFixedMerFreq/query.maxFixedMerFreq;
			int isGiantRepeat = ((queryType >> 1) & (targetType >> 1)) + 1;
			target.isHitchhiked = target.isHitchhiked || (query.isRepeat && freqDiff < (m_params.hhRatio/isGiantRepeat));	//HIGH --> LOW
			query.isHitchhiked = query.isHitchhiked || (target.isRepeat && freqDiff > (isGiantRepeat/m_params.hhRatio));	//LOW  --> HIGH
		}
	}
	SeedVector finalSeedVec, outcastSeedVec;
	finalSeedVec.reserve(initSeedVec.size());
	outcastSeedVec.reserve(initSeedVec.size() >> 1);
	for(const auto& iter : initSeedVec)
	{
		if(iter.isHitchhiked)
			outcastSeedVec.push_back(iter);
		else
			finalSeedVec.push_back(iter);
	}
	if(m_params.DebugSeed)
	{
		std::ostream* pOutcastSeedWriter = createWriter(m_params.directory + "seed/error/" + result.readid + ".seed");
		write(*pOutcastSeedWriter, outcastSeedVec);
		delete pOutcastSeedWriter;
	}
	return finalSeedVec;
}

void PacBioSelfCorrectionProcess::write(std::ostream& outfile, const SeedVector& seedVec) const
{
	for(const auto& iter : seedVec)
		outfile
		<< iter.seedStr << "\t"
		<< iter.maxFixedMerFreq << "\t"
		<< iter.seedStartPos << "\t"
		<< (iter.isRepeat ? "Yes" : "No") << "\n";
}
//Correct sequence by FMWalk & MSAlignment. Noted by KuanWeiLee
void PacBioSelfCorrectionProcess::initCorrect(std::string& readSeq, const SeedVector& seedVec, SeedVector& pieceVec, PacBioSelfCorrectionResult& result)
{
	if(m_params.OnlySeed || seedVec.size() < 2) return;
	pieceVec.push_back(seedVec[0]);
	pieceVec.back().seedStr.reserve(readSeq.length());
	std::ostream* pExtWriter = nullptr;
	std::ostream* pDpWriter = nullptr;
	std::ostream* pExtDebugFile = nullptr;
	if(m_params.DebugSeed)
	{
		pExtWriter    = createWriter(m_params.directory + "extend/" + result.readid + ".ext");
		pDpWriter     = createWriter(m_params.directory + "extend/" + result.readid + ".dp");
		pExtDebugFile = createWriter(m_params.directory + "extend/" + result.readid + ".fa");
	}

	int case_number = 1;
	for(SeedVector::const_iterator iterTarget = seedVec.begin() + 1; iterTarget != seedVec.end(); iterTarget++, case_number++)
	{
		int isFMExtensionSuccess = 0, firstFMExtensionType = 0;
		SeedFeature& source = pieceVec.back();
		std::string mergedSeq;

		for(int next = 0; next < m_params.numOfNextTarget && (iterTarget + next) != seedVec.end() ; next++)
		{
			const SeedFeature& target = *(iterTarget + next);
/*
			debugExtInfo debug(
								m_params.DebugExtend, pExtDebugFile, result.readid, case_number,
								source.seedEndPos - source.seedLength+1, source.seedEndPos,
								(*iterTarget).seedStartPos,(*iterTarget).seedEndPos, true
							);
			isFMExtensionSuccess = correctByFMExtension(source, target, readSeq, mergedSeq, result, debug);
/*/
			debugExtInfo debug;
			isFMExtensionSuccess = correctByFMExtension(source, target, readSeq, mergedSeq, result, debug);
//*/
			firstFMExtensionType = (next == 0 ? isFMExtensionSuccess : firstFMExtensionType);
			if(isFMExtensionSuccess > 0)
			{
				result.totalWalkNum++;
				source.append(mergedSeq, target);
				iterTarget += next;
				case_number+= next;
				break;
			}

		}

		if(isFMExtensionSuccess <= 0)
		{
			const SeedFeature& target = *iterTarget;
			switch(firstFMExtensionType)
			{
				case -1:
					result.highErrorNum++;
					break;
				case -2:
					result.exceedDepthNum++;
					break;
				case -3:
					result.exceedLeaveNum++;
					break;
				case -4:
					result.highErrorNum++;
					firstFMExtensionType = -1;
					break;
				default:
					std::cerr << "Does it really happen?\n";
					exit(EXIT_FAILURE);
			}
			if(m_params.DebugSeed)
				*pExtWriter << source.seedStartPos << "\t" << target.seedStartPos << "\t" << (firstFMExtensionType + 4) << "\n";
			result.totalWalkNum++;
			bool isMSAlignmentSuccess = correctByMSAlignment(source, target, readSeq, mergedSeq, result);
			if(isMSAlignmentSuccess)
				source.append(mergedSeq, target);
			else
			{
				if(m_params.DebugSeed)
					*pDpWriter << source.seedStartPos << "\t" << target.seedStartPos << "\n";
				if(m_params.Split)
					pieceVec.push_back(target);
				else
				{
					mergedSeq = readSeq.substr(source.seedEndPos + 1, target.seedEndPos - source.seedEndPos);
					source.append(mergedSeq, target);
				}
				result.correctedLen += target.seedStr.length();
			}
		}
	}

	delete pExtWriter;
	delete pDpWriter;
	delete pExtDebugFile;
}

int PacBioSelfCorrectionProcess::correctByFMExtension
(const SeedFeature& source, const SeedFeature& target, const std::string& in, std::string& out, PacBioSelfCorrectionResult& result, debugExtInfo& debug)
{
	int interval = target.seedStartPos - source.seedEndPos - 1;
	int extendKmerSize = MIN(source.endBestKmerSize, target.startBestKmerSize) - 2;
	if(source.isRepeat || target.isRepeat)
	{
		extendKmerSize = MIN(source.seedLength, target.seedLength);
		extendKmerSize = MIN(extendKmerSize, m_params.startKmerLength + 2);
	}
//	if(source.isRepeat && target.isRepeat)
//		extendKmerSize = MIN(source.seedLength, target.seedLength);
//	size_t extendKmerSize = m_params.startKmerLength;

	std::string src, trg, path;
	src = source.seedStr.substr(source.seedLength - extendKmerSize);
		debug.sourceReduceSize(source.seedLength - extendKmerSize);
	trg = target.seedStr;
	path = in.substr(source.seedEndPos + 1, interval);
	int min_SA_threshold = 3, isFMExtensionSuccess = 0;
	min_SA_threshold = m_params.PBcoverage > 60 ? ((m_params.PBcoverage / 60) * 3) : min_SA_threshold;
	bool isFromRtoU = source.isRepeat && !target.isRepeat;
	if(isFromRtoU)
	{
		std::swap(src,trg);
		src = reverseComplement(src);
		trg = reverseComplement(trg);
		path = reverseComplement(path);
			debug.reverseStrand();
	}

	Timer* FMTimer = new Timer("FM Time",true);
	FMWalkResult2 fmwalkresult;
	LongReadSelfCorrectByOverlap OverlapTree
	(src, path, trg, interval, extendKmerSize, extendKmerSize + 2, m_params.FM_params, min_SA_threshold, debug);
	isFMExtensionSuccess = OverlapTree.extendOverlap(fmwalkresult);
	result.Timer_FM += FMTimer->getElapsedWallTime();
	delete FMTimer;

	if(isFMExtensionSuccess < 0) return isFMExtensionSuccess;
	if(isFromRtoU)
	{
		fmwalkresult.mergedSeq = reverseComplement(fmwalkresult.mergedSeq);
		fmwalkresult.mergedSeq += reverseComplement(src).substr(extendKmerSize);
	}
	out = fmwalkresult.mergedSeq;
	out.erase(0,extendKmerSize);
	result.correctedLen += out.length();
	result.seedDis += interval;
	result.FMNum++;
	return isFMExtensionSuccess;
}

bool PacBioSelfCorrectionProcess::correctByMSAlignment
(const SeedFeature& source, const SeedFeature& target, const std::string& in, std::string& out, PacBioSelfCorrectionResult& result)
{
	if(m_params.NoDp) return false;
	int interval = target.seedStartPos - source.seedEndPos - 1;
	int extendKmerSize = MIN(source.endBestKmerSize, target.startBestKmerSize) - 2;
	if(source.isRepeat || target.isRepeat)
	{
		extendKmerSize = MIN(source.seedLength, target.seedLength);
		extendKmerSize = MIN(extendKmerSize, m_params.startKmerLength + 2);
	}
	std::string src, trg, path;
	src = source.seedStr.substr(source.seedLength - extendKmerSize);
	trg = target.seedStr;
	path = in.substr(source.seedEndPos + 1, interval);
	path = src + path + trg;
	double identity = 0.65;
	size_t totalMaxFixedMerFreq = source.maxFixedMerFreq + target.maxFixedMerFreq, min_call_coverage = 15;
	identity += (totalMaxFixedMerFreq > 50  ? 0.05 : 0);
	identity += (totalMaxFixedMerFreq > 100 ? 0.05 : 0);
	min_call_coverage = totalMaxFixedMerFreq > 50 ? totalMaxFixedMerFreq * 0.4 : min_call_coverage;

	Timer* DPTimer = new Timer("DP Time", true);
	MultipleAlignment maquery =
	LongReadOverlap::buildMultipleAlignment
	(path, extendKmerSize, extendKmerSize, path.length()/10, identity, m_params.PBcoverage, m_params.indices);
	result.Timer_DP += DPTimer->getElapsedWallTime();
	delete DPTimer;


	if(maquery.getNumRows() <= 3) return false;
	out = maquery.calculateBaseConsensus(min_call_coverage, -1);
	out.erase(0,extendKmerSize);
	result.correctedLen += out.length();
	result.seedDis += interval;
	result.DPNum++;
	return true;
}

//
//
//
PacBioSelfCorrectionPostProcess::PacBioSelfCorrectionPostProcess(
		std::string correctFile,
		std::string discardFile,
		const PacBioSelfCorrectionParameters params)
:	m_params(params),
	m_totalReadsLen(0),
	m_correctedLen(0),
	m_totalSeedNum(0),
	m_totalWalkNum(0),
	m_highErrorNum(0),
	m_exceedDepthNum(0),
	m_exceedLeaveNum(0),
	m_FMNum(0),
	m_DPNum(0),
	m_OutcastNum(0),
	m_seedDis(0),
	m_Timer_Seed(0),
	m_Timer_FM(0),
	m_Timer_DP(0)
{
	m_pCorrectWriter = createWriter(correctFile);
	m_pDiscardWriter = createWriter(discardFile);
}

//
PacBioSelfCorrectionPostProcess::~PacBioSelfCorrectionPostProcess()
{
	m_OutcastNum = m_totalWalkNum - m_FMNum - m_DPNum;
	if(m_totalWalkNum>0 && m_totalReadsLen>0)
	{
		std::cout << "\n"
		<< "TotalReadsLen: " << m_totalReadsLen << "\n"
		<< "CorrectedLen: " << m_correctedLen << ", ratio: " << (float)(m_correctedLen)/m_totalReadsLen << "\n"
		<< "TotalSeedNum: " << m_totalSeedNum << "\n"
		<< "TotalWalkNum: " << m_totalWalkNum << "\n"
		<< "FMNum: " << m_FMNum << ", ratio: " << (float)(m_FMNum * 100)/m_totalWalkNum << "%\n"
		<< "DPNum: " << m_DPNum << ", ratio: " << (float)(m_DPNum * 100)/m_totalWalkNum << "%\n"
        << "OutcastNum: " << m_OutcastNum << ", ratio: " << (float)(m_OutcastNum * 100)/m_totalWalkNum << "%\n"
		<< "HighErrorNum: " << m_highErrorNum << ", ratio: " << (float)(m_highErrorNum * 100)/(m_DPNum + m_OutcastNum) << "%\n"
		<< "ExceedDepthNum: " << m_exceedDepthNum << ", ratio: " << (float)(m_exceedDepthNum * 100)/(m_DPNum + m_OutcastNum) << "%\n"
		<< "ExceedLeaveNum: " << m_exceedLeaveNum << ", ratio: " << (float)(m_exceedLeaveNum * 100)/(m_DPNum + m_OutcastNum) << "%\n"
		<< "DisBetweenSeeds: " << m_seedDis/m_totalWalkNum << "\n"
        << "Time of searching Seeds: " << m_Timer_Seed  << "\n"
        << "Time of searching FM: " << m_Timer_FM  << "\n"
        << "Time of searching DP: " << m_Timer_DP  << "\n";
	}
	delete m_pCorrectWriter;
	delete m_pDiscardWriter;
}


// Writting results for kmerize and validate
void PacBioSelfCorrectionPostProcess::process(const SequenceWorkItem& item, const PacBioSelfCorrectionResult& result)
{
	if(result.merge)
	{
		m_totalReadsLen += result.totalReadsLen;
		m_correctedLen += result.correctedLen;
		m_totalSeedNum += result.totalSeedNum;
		m_totalWalkNum += result.totalWalkNum;
		m_highErrorNum += result.highErrorNum;
		m_exceedDepthNum += result.exceedDepthNum;
		m_exceedLeaveNum += result.exceedLeaveNum;
		m_FMNum += result.FMNum;
        m_DPNum += result.DPNum;
		m_seedDis += result.seedDis;
        m_Timer_Seed += result.Timer_Seed;
        m_Timer_FM += result.Timer_FM;
        m_Timer_DP += result.Timer_DP;

		for(std::vector<DNAString>::const_iterator iter = result.correctedStrs.begin(); iter != result.correctedStrs.end(); iter++)
		{
			size_t i = iter - result.correctedStrs.begin();
			SeqItem mergeRecord;
			std::stringstream ss;
			ss << item.read.id << "_" << i << (*iter).toString().length();
			mergeRecord.id = ss.str();
			mergeRecord.seq = *iter;
			mergeRecord.write(*m_pCorrectWriter);
		}
	}
	else
	{
		// write into discard.fa
		SeqItem mergeRecord;
		mergeRecord.id = item.read.id;
		mergeRecord.seq = item.read.seq;
		mergeRecord.write(*m_pDiscardWriter);
	}
}