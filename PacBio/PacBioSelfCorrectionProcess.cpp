///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioSelfCorrectionProcess.cpp - Self-correction using FM-index walk for PacBio reads
//
#include <algorithm>
#include <memory>
#include "PacBioSelfCorrectionProcess.h"
#include "LongReadProbe.h"
#include "LongReadOverlap.h"
#include "Util.h"
#include "Timer.h"
#include "KmerFeature.h"

// PacBio Self Correction by Ya and YTH, v20151202.
// 1. Identify highly-accurate seeds within PacBio reads
// 2. For each pair of seeds, perform kmer extension using local kmer frequency collected by FM-index extension
PacBioSelfCorrectionResult PacBioSelfCorrectionProcess::process(const SequenceWorkItem& workItem)
{
	PacBioSelfCorrectionResult result;
    result.readid = workItem.read.id;
	std::string readSeq = workItem.read.seq.toString();
	const size_t readSeqLen = readSeq.length();
	SeedFeature::SeedVector seedVec, pieceVec;

	//allocate space for kmers on each position
	for(auto& iter : m_params.kmerPool)
		KmerFeature::kmerRec[iter] = std::unique_ptr<KmerFeature[]>(new KmerFeature[readSeqLen]);

	//Part 1: start searching seeds
    Timer* seedTimer = new Timer("Seed Time", true);
	LongReadProbe::readid = result.readid;
	LongReadProbe::searchSeedsWithHybridKmers(readSeq, seedVec);
	result.totalSeedNum = seedVec.size();
	result.Timer_Seed = seedTimer->getElapsedWallTime();
	delete seedTimer;

	//Part 2:start correcting sequence
    initCorrect(readSeq, seedVec, pieceVec, result);

	//free space for kmers on each position
	KmerFeature::kmerRec.clear();

	result.merge = !pieceVec.empty();
	result.totalReadsLen = readSeq.length();
	for(const auto& iter : pieceVec)
		result.correctedStrs.push_back(iter.seedStr);
	return result;
}
//Correct sequence by FMWalk & MSAlignment; it's a workflow control module. Noted by KuanWeiLee 18/3/12
void PacBioSelfCorrectionProcess::initCorrect(std::string& readSeq, const SeedFeature::SeedVector& seedVec, SeedFeature::SeedVector& pieceVec, PacBioSelfCorrectionResult& result)
{
	std::ostream* pExtWriter = nullptr;
	std::ostream* pDpWriter = nullptr;
	std::ostream* pExtDebugFile = nullptr;
	if(m_params.OnlySeed || seedVec.size() < 2) return;

	//push first seed into vector and reserve space for fast expansion
	pieceVec.push_back(seedVec[0]);
	pieceVec.back().seedStr.reserve(readSeq.length());
	if(m_params.DebugSeed)
	{
		pExtWriter    = createWriter(m_params.directory + "extend/" + result.readid + ".ext");
		pDpWriter     = createWriter(m_params.directory + "extend/" + result.readid + ".dp");
	}
/*
	if(m_params.DebugExtend)
		pExtDebugFile = createWriter(m_params.directory + "extensionFile/" + result.readid + ".fa");
//*/
	int case_number = 1;
	for(SeedFeature::SeedVector::const_iterator iterTarget = seedVec.begin() + 1; iterTarget != seedVec.end(); iterTarget++, case_number++)
	{
		int isFMExtensionSuccess = 0, firstFMExtensionType = 0;
		SeedFeature& source = pieceVec.back();
		std::string mergedSeq;

		for(int next = 0; next < m_params.numOfNextTarget && (iterTarget + next) != seedVec.end() ; next++)
		{
			const SeedFeature& target = *(iterTarget + next);
/*
			debugExtInfo debug( m_params.DebugExtend, pExtDebugFile, result.readid, case_number);
/*/
			debugExtInfo debug;
//*/

			isFMExtensionSuccess = correctByFMExtension(source, target, readSeq, mergedSeq, result, debug);
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
					mergedSeq = readSeq.substr((source.seedEndPos + 1), (target.seedEndPos - source.seedEndPos));
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
	bool isFromRtoU = source.isRepeat && !target.isRepeat;
	int isFMExtensionSuccess = 0;

	// Set initial extend size
		int initExtSize = std::min(source.endBestKmerSize, target.startBestKmerSize) - 2;
			if(source.isRepeat || target.isRepeat)
			{
				initExtSize = std::min(source.seedLen, target.seedLen);
				initExtSize = std::min(initExtSize, m_params.startKmerLen + 2);
			}

	// Set seed Information
		seedPair extSeeds(source, target);
		extSeeds.reduceSourceBy(source.seedLen - initExtSize);

	// Set the raw read path which shall be corrected.
		int seedDistance     = target.seedStartPos - source.seedEndPos - 1;
		std::string readPath = in.substr(source.seedEndPos + 1, seedDistance);

		if(isFromRtoU)
		{
			extSeeds.reverseSeed();
			readPath = reverseComplement(readPath);
		}

	// Set the normal threshold by extension
		int min_SA_threshold;
		if(m_params.PBcoverage > 60)
			min_SA_threshold = (m_params.PBcoverage / 60) * 3;
		else
			min_SA_threshold = 3;

	// FM extension
		Timer* FMTimer = new Timer("FM Time",true);
			FMWalkResult2 fmwalkresult;
			LongReadSelfCorrectByOverlap OverlapTree
				(extSeeds, readPath, seedDistance, initExtSize, initExtSize + 2, m_params.FM_params, min_SA_threshold, debug);
			isFMExtensionSuccess = OverlapTree.extendOverlap(fmwalkresult);
			result.Timer_FM += FMTimer->getElapsedWallTime();
		delete FMTimer;

	if(isFMExtensionSuccess < 0)
		return isFMExtensionSuccess;

	// Recover reverse Sequence
		if(isFromRtoU)
		{
			fmwalkresult.mergedSeq  = reverseComplement(fmwalkresult.mergedSeq);
			fmwalkresult.mergedSeq += reverseComplement(extSeeds.source.seq).substr(initExtSize);
		}

	// Set the output sequence
		out = fmwalkresult.mergedSeq;
		out.erase(0,initExtSize);
		result.correctedLen += out.length();
		result.seedDis      += seedDistance;
		result.FMNum++;

	return isFMExtensionSuccess;

}

bool PacBioSelfCorrectionProcess::correctByMSAlignment
(const SeedFeature& source, const SeedFeature& target, const std::string& in, std::string& out, PacBioSelfCorrectionResult& result)
{
	if(m_params.NoDp) return false;
	int interval = target.seedStartPos - source.seedEndPos - 1;
	int extendKmerSize = std::min(source.endBestKmerSize, target.startBestKmerSize) - 2;
	if(source.isRepeat || target.isRepeat)
	{
		extendKmerSize = std::min(source.seedLen, target.seedLen);
		extendKmerSize = std::min(extendKmerSize, m_params.startKmerLen + 2);
	}
	std::string src, trg, path;
	src = source.seedStr.substr(source.seedLen - extendKmerSize);
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
		const PacBioSelfCorrectionParameters& params)
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
		<< "FMNum: "        << m_FMNum      << ", ratio: " << (float)(m_FMNum      * 100)/m_totalWalkNum << "%\n"
		<< "DPNum: "        << m_DPNum      << ", ratio: " << (float)(m_DPNum      * 100)/m_totalWalkNum << "%\n"
        << "OutcastNum: "   << m_OutcastNum << ", ratio: " << (float)(m_OutcastNum * 100)/m_totalWalkNum << "%\n"
		<< "HighErrorNum: "    << m_highErrorNum   << ", ratio: " << (float)(m_highErrorNum   * 100)/(m_DPNum + m_OutcastNum) << "%\n"
		<< "ExceedDepthNum: "  << m_exceedDepthNum << ", ratio: " << (float)(m_exceedDepthNum * 100)/(m_DPNum + m_OutcastNum) << "%\n"
		<< "ExceedLeaveNum: "  << m_exceedLeaveNum << ", ratio: " << (float)(m_exceedLeaveNum * 100)/(m_DPNum + m_OutcastNum) << "%\n"
		<< "DisBetweenSeeds: " << m_seedDis/m_totalWalkNum << "\n"
        << "Time of searching Seeds: " << m_Timer_Seed  << "\n"
        << "Time of searching FM: "    << m_Timer_FM  << "\n"
        << "Time of searching DP: "    << m_Timer_DP  << "\n";
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
			size_t index = iter - result.correctedStrs.begin();
			SeqItem mergeSeq;
			std::string flag = m_params.Split ? ("_" + std::to_string(index)) : "";
			mergeSeq.id = item.read.id + flag;
			mergeSeq.seq = *iter;
			mergeSeq.write(*m_pCorrectWriter);
		}
	}
	else
	{
		// write into discard.fa
		SeqItem mergeSeq;
		mergeSeq.id = item.read.id;
		mergeSeq.seq = item.read.seq;
		mergeSeq.write(*m_pDiscardWriter);
	}
}