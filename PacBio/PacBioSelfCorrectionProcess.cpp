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
#include <iomanip>
#include <time.h>
#include "SAIPBHybridCTree.h"
#include "LongReadOverlap.h"
#include "Timer.h"
#include "KmerThresholdTable.h"

#define OVERTHRESHOLD( a,b,c,d ) ((a >= b) && (c >= 1) && (d >= 1))
#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)

// PacBio Self Correction by Ya and YTH, v20151202.
// 1. Identify highly-accurate seeds within PacBio reads
// 2. For each pair of seeds, perform kmer extension using local kmer frequency collected by FM-index extension
PacBioSelfCorrectionResult PacBioSelfCorrectionProcess::process(const SequenceWorkItem& workItem)
{	
	PacBioSelfCorrectionResult result;
    result.readid = workItem.read.id;
	result.kd.setID(result.readid);
	
	std::string readSeq = workItem.read.seq.toString();	
	SeedVector seedVec, pieceVec;
	searchSeedsWithHybridKmers(readSeq, seedVec, result);
	//Push the first seed into pieceVec, which will be popped later as source seed

	//Give up reads with less than 2 seeds
	if(m_params.OnlySeed || seedVec.size() < 2) return result;
	result.correctedLen += seedVec[0].seedStr.length();
	pieceVec.push_back(seedVec[0]);
    //Reserve sufficient str length for fast append
    pieceVec.back().seedStr.reserve(readSeq.length());
    initCorrect(readSeq, seedVec, pieceVec, result);
	
	result.merge = true;
	result.totalReadsLen = readSeq.length();
	for(const auto& iter : pieceVec)
		result.correctedStrs.push_back(iter.seedStr);
	return result;
    
}

// Search seeds with static and dynamic kmers. Noted by KuanWeiLee 20171027
void PacBioSelfCorrectionProcess::searchSeedsWithHybridKmers(const std::string& readSeq, SeedVector& seedVec, PacBioSelfCorrectionResult& result)
{
    Timer* seedTimer = new Timer("Seed Time",true);
	int type[readSeq.length()]{0};
	unsigned int offset = determineTableType(readSeq, type, result);
    const unsigned int staticKmerSize = m_params.kmerLength - offset;
	if(readSeq.length() < staticKmerSize) return;
    
	//Part 1 : Slide the sequence with a static-kmer. Noted by KuanWeiLee
	BistrandBWTInterval FixedMerInterval[readSeq.length() - staticKmerSize + 1];
	//std::ostream* pKfWriter = createWriter(m_params.directory + "seed/." + result.readid + ".kf");
    for(size_t i = 0 ; i <= (readSeq.length() - staticKmerSize) ; i++)
	{
        std::string kmer = readSeq.substr(i, staticKmerSize);
		BistrandBWTInterval object;
		object.fwdInterval = BWTAlgorithms::findInterval(m_params.indices.pRBWT, reverse(kmer));
		object.rvcInterval = BWTAlgorithms::findInterval(m_params.indices.pBWT, reverseComplement(kmer));
		FixedMerInterval[i] = object;
		int kmerFreq = object.getFreq();
		result.kd.add(kmerFreq);
	//	*pKfWriter << i + 1 << "\t" << kmerFreq << "\n";
    }
	result.kd.computeKDAttributes(1);
	//delete pKfWriter;
	
	float** table = KmerThresholdTable::m_table;
	//Part 2 : Search seeds; slide through the read sequence with hybrid-kmers. Noted by KuanWeiLee
	//kmer[Start/Move]Pos indicate the starting & moving positions of the static-kmer.
	for(size_t kmerStartPos = 0; kmerStartPos <= (readSeq.length() - staticKmerSize); kmerStartPos++)
	{
		bool isSeed = false, isRepeat = false;
		std::string kmer = readSeq.substr(kmerStartPos,staticKmerSize);
		unsigned int dynamicKmerSize = staticKmerSize;
		BWTInterval fwdInterval = FixedMerInterval[kmerStartPos].fwdInterval;
		BWTInterval rvcInterval = FixedMerInterval[kmerStartPos].rvcInterval;
		size_t fwdKmerFreq, rvcKmerFreq, dynamicKmerFreq;
		size_t staticKmerFreq, maxFixedMerFreq = FixedMerInterval[kmerStartPos].getFreq();
		float staticKmerThresholdValue, dynamicKmerThresholdValue, repeatKmerThresholdValue;
		size_t seedStartPos = kmerStartPos, seedEndPos;
		float GCratio, freqsDiff;
		for(size_t kmerMovePos = kmerStartPos; kmerMovePos <= (readSeq.length() - staticKmerSize); kmerMovePos++)
		{	
			
			if(isSeed)
			{
				char s = readSeq[kmerMovePos + staticKmerSize - 1];//s:step
				char rcs = complement(s);
				kmer += s;
				dynamicKmerSize++;
				BWTAlgorithms::updateInterval(fwdInterval, s,m_params.indices.pRBWT);
				BWTAlgorithms::updateInterval(rvcInterval, rcs,m_params.indices.pBWT);
			}
			staticKmerThresholdValue = table[type[kmerMovePos]][staticKmerSize];
			dynamicKmerThresholdValue = table[type[seedStartPos]][dynamicKmerSize];
			repeatKmerThresholdValue = type[kmerMovePos] == 2 ? staticKmerThresholdValue : staticKmerThresholdValue * 5;
			fwdKmerFreq = fwdInterval.getFreq();
			rvcKmerFreq = rvcInterval.getFreq();			
			dynamicKmerFreq = fwdKmerFreq + rvcKmerFreq;
			staticKmerFreq = FixedMerInterval[kmerMovePos].getFreq();
			freqsDiff = (float)staticKmerFreq/maxFixedMerFreq;
			//The order between general seed extentsion and kmer hitchhike strategy may make difference 
			//,which needs furher observation. Noted by KuanWeiLee 20171027
			//Gerneral seed extension strategy.
			if  (
				   dynamicKmerSize > m_params.kmerLengthUpperBound											//1.over length
				|| staticKmerFreq < staticKmerThresholdValue												//2.fixed kmer frequency
				|| !OVERTHRESHOLD(dynamicKmerFreq, dynamicKmerThresholdValue, fwdKmerFreq, rvcKmerFreq)		//3.dynamic kmer frequency
				)
			{
				if(isSeed) kmer.erase(--dynamicKmerSize);
				break;
			}
			//Kmer Hitchhike strategy.
			if(isRepeat && freqsDiff < m_params.khhRatio)													//4.hitchhiking kmer (HIGH-->LOW)
			{
				kmer.erase(--dynamicKmerSize);
				kmerStartPos++;
				break;
			}
			else if(staticKmerFreq >= repeatKmerThresholdValue && freqsDiff > (1/m_params.khhRatio))		//4.hitchhiking kmer (LOW-->HIGH)
			{
				isSeed = false;
				kmerStartPos = kmerMovePos - 1;
				break;
			}
			isSeed = true;
			seedEndPos = seedStartPos + dynamicKmerSize - 1;
			kmerStartPos = seedEndPos;
			isRepeat = isRepeat || (staticKmerFreq >= repeatKmerThresholdValue);
			maxFixedMerFreq = MAX(maxFixedMerFreq, staticKmerFreq);
		}
		if(isSeed)
		{
			//Low Complexity strategy.
			if(isLowComplexity(kmer, GCratio)) continue;
			SeedFeature newSeed(kmer, seedStartPos, maxFixedMerFreq, isRepeat, staticKmerSize, m_params.PBcoverage);
			newSeed.estimateBestKmerSize(m_params.indices);
			seedVec.push_back(newSeed);
		}
	}
	//Seed Hitchhike strategy.
	seedVec = removeHitchhikingSeeds(seedVec, type, result);
	//[Debugseed] Output seeds & kd statistics for each reads.Noted by KuanWeiLee
    if(m_params.DebugSeed)
    {
        std::ostream* pSeedWriter = createWriter(m_params.directory + "seed/" + result.readid + ".seed");
		write(*pSeedWriter, seedVec);
		delete pSeedWriter;
		
		std::ostream* pKdWriter = createWriter(m_params.directory + "kmer-stat", std::ios_base::app);
		result.kd.write(*pKdWriter,KmerDistribution::TYPE::ATTRIBUTE);
		delete pKdWriter;
	}
	
	result.totalSeedNum = seedVec.size();
	result.Timer_Seed = seedTimer->getElapsedWallTime(); 
    delete seedTimer;
}
//KmerThresholdTable type is set dynamically using a sliding fixed-mer on each position of the sequence.
//Noted by KuanWeiLee 20180118
unsigned int PacBioSelfCorrectionProcess::determineTableType(const std::string& seq, int* const type, const PacBioSelfCorrectionResult& result)
{
	std::ostream* totalWriter = createWriter(m_params.directory + "total", std::ios_base::app);
	//std::ostream* recordWriter = createWriter(m_params.directory + "seed/." + result.readid + ".rec");
	//std::ostream* ratioWriter = createWriter(m_params.directory + "seed/." + result.readid + ".rr");
	unsigned int offset = 0;
	size_t seqlen = seq.length();
	std::fill_n(type, seqlen, 1);
	int record[seqlen] = {0};
	for(size_t i = 0; i <= (seqlen - m_params.scanningKmerLength); i++)
	{
		std::string kmer = seq.substr(i, m_params.scanningKmerLength);
		int kmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
		record[i] = kmerFreq;
		//*recordWriter << i << "\t" << kmerFreq << "\n";
	}
	int range = 300;
	int x = m_params.PBcoverage;
	int y = m_params.scanningKmerLength;
//	float lowcov = KmerThresholdTable::calculate(0, x, y);
//	float unique = KmerThresholdTable::calculate(1, x, y);
	float repeat = KmerThresholdTable::calculate(2, x ,y);
	int head = 0, tail = -1;
	int leftmost = seqlen - 1, rightmost = 0, repeatcount = 0;
	std::map<std::string, int> set;
	for(size_t i = 0; i < seqlen; i++)
	{
		int left = i - (range >> 1);
		left = MAX(left, 0);
		int right = i + (range >> 1);
		right = MIN(right, (int)(seqlen - 1));
		while(tail < right)
		{
			int freq = record[++tail];
			if(freq == 0)
			{
				set["empty"]++;
				continue;
			}
			if(freq >= repeat)
			{
				set["repeat"]++;
				continue;
			}
			set["other"]++;
		}
		while(head < left)
		{
			int freq = record[head++];
			if(freq >= repeat)
			{
				set["repeat"]--;
				continue;
			}
			set["other"]--;
		}
		int size = (right - left + 1) - set["empty"];
		float ratio = (float)set["repeat"]/size + 0.0005;
		if(ratio >= 0.02)
		{
			type[i] = 2;
			repeatcount++;
			leftmost = MIN(leftmost, (int)i);
			rightmost = MAX(rightmost, (int)i);
		}
		//*ratioWriter << i << "\t" << ratio << "\n";
	}
	*totalWriter << result.readid << "\t" << (float)leftmost/seqlen << "\t" << (float)(seqlen - rightmost)/seqlen << "\t" << (float)repeatcount/seqlen << "\n";
	if	(
		   (float)repeatcount/seqlen >= 0.5
		&& (float)(leftmost + (seqlen -rightmost))/seqlen <= 0.1
		)
	{
		std::fill_n(type, seqlen, 2);
		offset = 4;
	}
	//delete totalWriter;
	//delete recordWriter;
	//delete ratioWriter;
	return offset;
}
//Kmer & Seed Hitchhike strategy would maitain seed-correctness, 
//once the sequence is stuck between the ambiguity from uniqu to repeat mode.
//Noted by KuanWeiLee 20180106
PacBioSelfCorrectionProcess::SeedVector PacBioSelfCorrectionProcess::removeHitchhikingSeeds(SeedVector initSeedVec, int const *type, PacBioSelfCorrectionResult& result)
{
	int x = m_params.PBcoverage;
	int y = m_params.kmerLength;
	float overfrequency = KmerThresholdTable::calculate(2, x ,y) * 5;
	for(SeedVector::iterator iterQuery = initSeedVec.begin(); iterQuery != initSeedVec.end(); iterQuery++)
	{
		SeedFeature& query = *iterQuery;
		SeedVector::iterator iterTarget = iterQuery + 1;
		if(type[query.seedStartPos] == 2 && query.maxFixedMerFreq >= overfrequency) continue;
		//if(query.isHitchhiked) continue;
		while(iterTarget != initSeedVec.end() && ((*iterTarget).seedStartPos - query.seedEndPos) <= m_params.repaetDistance)
		{
			SeedFeature& target = *iterTarget;
			iterTarget++;
			if(type[target.seedStartPos] == 2 && target.maxFixedMerFreq >= overfrequency) continue;
			//if(target.isHitchhiked) continue;
			float freqsDiff = (float)target.maxFixedMerFreq/query.maxFixedMerFreq;
			int isGiantRepeat = (type[query.seedStartPos] ==2 && type[target.seedStartPos] == 2) ? 2 : 1;
			target.isHitchhiked = target.isHitchhiked || (query.isRepeat && freqsDiff < (m_params.shhRatio/isGiantRepeat));	//HIGH --> LOW
			query.isHitchhiked = query.isHitchhiked || (target.isRepeat && freqsDiff > (isGiantRepeat/m_params.shhRatio));	//LOW  --> HIGH
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
		std::ostream* pOutcastSeedWriter = createWriter(m_params.directory + "seed/shh/" + result.readid + ".seed");
		write(*pOutcastSeedWriter, outcastSeedVec);
		delete pOutcastSeedWriter;
	}
	return finalSeedVec;
}

bool PacBioSelfCorrectionProcess::isLowComplexity (const std::string& seq, float& GCratio, float threshold)
{
	size_t seqlen = seq.length();
	int count[4]{0};
	for (size_t i = 0; i < seqlen; i++)
	{
		switch(seq[i])
		{
			case 'A': count[0]++; break;
			case 'T': count[1]++; break;
			case 'C': count[2]++; break;
			case 'G': count[3]++; break;
			default: assert(false);
		}
	}
	GCratio = (float)(count[2] + count[3])/seqlen;
	bool case1 = false, case2 = false;
	int num = 0;
	for(const auto& iter : count)
	{
		case1 = case1 || ((float)iter/seqlen >= threshold);
		num += (iter == 0 ? 1 : 0);
	}
	case2 = (num == 2);
	return case1 || case2;
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
	std::ostream* pExtendWriter = nullptr;
	if(m_params.DebugExtend)  pExtendWriter = createWriter(m_params.directory + "extend/" + result.readid + ".ext");
	for(SeedVector::const_iterator iterTarget = seedVec.begin() + 1; iterTarget != seedVec.end(); iterTarget++)
	{
		int isFMExtensionSuccess = 0, firstFMExtensionType = 0;
		SeedFeature& source = pieceVec.back();
		std::string mergedSeq;
		
		for(int next = 0; next < m_params.numOfNextTarget && (iterTarget + next) != seedVec.end() ; next++)
		{
			const SeedFeature& target = *(iterTarget + next);
			isFMExtensionSuccess = correctByFMExtension(source, target, readSeq, mergedSeq, result);
			firstFMExtensionType = (next == 0 ? isFMExtensionSuccess : firstFMExtensionType);
			if(isFMExtensionSuccess > 0)
			{
				result.totalWalkNum++;
				source.append(mergedSeq, target);
				iterTarget += next;
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
					std::cout << "Does it really happen?\n";
					exit(EXIT_FAILURE);
			}
			if(m_params.DebugExtend)
				*pExtendWriter << source.seedStartPos << "\t" << target.seedStartPos << "\t" << (firstFMExtensionType + 4) << "\n";
			result.totalWalkNum++;
			bool isMSAlignmentSuccess = correctByMSAlignment(source, target, readSeq, mergedSeq, result);
			if(isMSAlignmentSuccess)
				source.append(mergedSeq, target);
			else if(!m_params.isSplit)
			{
				mergedSeq = readSeq.substr(source.seedEndPos + 1, target.seedEndPos - source.seedEndPos);
				source.append(mergedSeq, target);
				result.correctedLen += target.seedStr.length();
			}
			else
			{
				pieceVec.push_back(target);
				result.correctedLen += target.seedStr.length();
			}
		}
	}
	if(pExtendWriter != nullptr) delete pExtendWriter;
}

int PacBioSelfCorrectionProcess::correctByFMExtension
(const SeedFeature& source, const SeedFeature& target, const std::string& in, std::string& out, PacBioSelfCorrectionResult& result)
{
	int interval = target.seedStartPos - source.seedEndPos - 1;
	size_t extendKmerSize = MIN(source.endBestKmerSize, target.startBestKmerSize) - 2;
	if(source.isRepeat || target.isRepeat)
	{
		extendKmerSize = MIN(source.seedLength, target.seedLength);
		extendKmerSize = MIN(extendKmerSize, m_params.kmerLength + 2);
	}
	std::string src, trg, path;
	src = source.seedStr.substr(source.seedLength - extendKmerSize);
	trg = target.seedStr;
	path = in.substr(source.seedEndPos + 1, interval);
	int min_SA_threshold = 3, isFMExtensionSuccess = 0;
	min_SA_threshold = m_params.PBcoverage > 60 ? ((m_params.PBcoverage / 60) * 3) : min_SA_threshold;
	//float freqsDiff = (float)target.maxFixedMerFreq/source.maxFixedMerFreq;
	bool isFromRtoU = false;
	//isFromRtoU = source.isRepeat && freqsDiff < m_params.shhRatio;//HIGH --> LOW
	//isFromRtoU = source.isRepeat && source.maxFixedMerFreq > target.maxFixedMerFreq;
	isFromRtoU = source.isRepeat && !target.isRepeat;
	if(isFromRtoU)
	{
		std::swap(src,trg);
		src = reverseComplement(src);
		trg = reverseComplement(trg);
		path = reverseComplement(path);
	}

	Timer* FMTimer = new Timer("FM Time",true);	
	FMWalkResult2 fmwalkresult;
	LongReadSelfCorrectByOverlap OverlapTree
	(src, path, trg, interval, extendKmerSize, extendKmerSize + 2, m_params.FM_params, min_SA_threshold);
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
	size_t extendKmerSize = MIN(source.endBestKmerSize, target.startBestKmerSize) - 2;
	if(source.isRepeat || target.isRepeat)
	{
		extendKmerSize = MIN(source.seedLength, target.seedLength);
		extendKmerSize = MIN(extendKmerSize, m_params.kmerLength + 2);
	}
	std::string src, trg, path;
	src = source.seedStr.substr(source.seedLength - extendKmerSize);
	trg = target.seedStr;
	path = in.substr(source.seedEndPos + 1, interval);
	path = src + path + trg;
	double identity = 0.65;
	size_t totalMaxFixedMerFreq = source.maxFixedMerFreq + target.maxFixedMerFreq, min_call_coverage = 15;
	identity += (totalMaxFixedMerFreq > 50  ? 0.5 : 0);
	identity += (totalMaxFixedMerFreq > 100 ? 0.5 : 0);
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
	if(m_params.DebugSeed)
		m_pKdWriter = createWriter(m_params.directory + "kmer-dist");
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
	if(m_params.DebugSeed)
	{
		m_kd.write(*m_pKdWriter);
		delete m_pKdWriter;
	}
	delete m_pCorrectWriter;
	delete m_pDiscardWriter;
}


// Writting results for kmerize and validate
void PacBioSelfCorrectionPostProcess::process(const SequenceWorkItem& item, const PacBioSelfCorrectionResult& result)
{	
	if (result.merge)
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
	m_kd += result.kd;
}