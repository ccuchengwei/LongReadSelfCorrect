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

//Generic judgment of a kmer ( a,b,c,d ) = ( currentKmerFreqs,dynamicKmerThresholdValue,fwdKmerFreqs,rvcKmerFreqs )
#define OVERTHRESHOLD( a,b,c,d ) ( (a>=b) && (c>=1) && (d>=1) )
#define MIN(a,b) ( a<b ? a : b)
#define MAX(a,b) ( a>b ? a : b)

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
	if(seedVec.size() >= 2 && !m_params.OnlySeed)
	{
		result.correctedLen += seedVec[0].seedStr.length();		
		pieceVec.push_back(seedVec[0]);		
	}
	else
		//Give up reads with less than 2 seeds
		return result;
    //Reserve sufficient str length for fast append
    pieceVec.back().seedStr.reserve(readSeq.length());
    initCorrect(readSeq, seedVec, pieceVec, result);
	
	result.merge = true;
	result.totalReadsLen = readSeq.length();
	for(SeedVector::iterator iter = pieceVec.begin(); iter != pieceVec.end(); iter++)
		result.correctedStrs.push_back((*iter).seedStr);
	
	return result;
    
}

// Search seeds with fixed and dynamic kmer size. Noted by KuanWeiLee 20171027
void PacBioSelfCorrectionProcess::searchSeedsWithHybridKmers(const std::string& readSeq, SeedVector& seedVec, PacBioSelfCorrectionResult& result)
{
    Timer* seedTimer = new Timer("Seed Time",true);
    const unsigned int staticKmerSize = m_params.kmerLength;
	if(readSeq.length() < staticKmerSize) return;
    
	//Part 1 : Slide with a fixed kmer. Noted by KuanWeiLee
	std::vector<BistrandBWTInterval> FixedMerInterval;
	
    for(size_t i = 0 ; i <= readSeq.length() - staticKmerSize ; i++)
	{
        std::string kmer = readSeq.substr(i,staticKmerSize);
		BistrandBWTInterval object;
		object.fwdInterval = BWTAlgorithms::findInterval(m_params.indices.pRBWT, reverse(kmer));
		object.rvcInterval = BWTAlgorithms::findInterval(m_params.indices.pBWT, reverseComplement(kmer));
		FixedMerInterval.push_back(object);
		size_t currentKmerFreqs = object.getFreqs();
		result.kd.add(currentKmerFreqs);
    }
	result.kd.computeKDAttributes(1);
	//LOWCOV:1 UNIQUE:2 REPEAT:3
	//KmerThresholdTable::TYPE mode = KmerThresholdTable::TYPE::LOWCOV;
    KmerThresholdTable::TYPE mode = KmerThresholdTable::TYPE::UNIQUE;
	//KmerThresholdTable::TYPE mode = KmerThresholdTable::TYPE::REPEAT;
	
	float staticKmerThresholdValue, repeatKmerThresholdValue;
	staticKmerThresholdValue = KmerThresholdTable::get(staticKmerSize,mode);
	repeatKmerThresholdValue = mode == 3 ? staticKmerThresholdValue : staticKmerThresholdValue * 5;
	
	//Part 2 : Search seeds; slide through the read sequence with dynamic kmers. Noted by KuanWeiLee
	for(size_t startPos = 0; startPos <= readSeq.length() - staticKmerSize; startPos++)
	{
		bool isSeed = false, isRepeat = false;
		std::string kmer = readSeq.substr(startPos,staticKmerSize);
		unsigned int dynamicKmerSize = staticKmerSize;
		BWTInterval fwdInterval = FixedMerInterval[startPos].fwdInterval,
					rvcInterval = FixedMerInterval[startPos].rvcInterval;
		size_t fwdKmerFreqs, rvcKmerFreqs, currentKmerFreqs;
		size_t fixedMerFreqs, maxFixedMerFreqs = FixedMerInterval[startPos].getFreqs();
		float dynamicKmerThresholdValue;
		size_t seedStartPos = startPos, seedEndPos;
		float GCratio, freqsDiff;
		for(size_t movePos = startPos; movePos <= readSeq.length() - staticKmerSize; movePos++, isSeed = true)
		{	
			
			if(isSeed)
			{
				char s = readSeq[movePos + staticKmerSize - 1];//s:step
				char rcs = complement(s);
				kmer += s;
				dynamicKmerSize++;
				BWTAlgorithms::updateInterval(fwdInterval,s,m_params.indices.pRBWT);
				BWTAlgorithms::updateInterval(rvcInterval,rcs,m_params.indices.pBWT);
			}
			isRepeat = isRepeat || (maxFixedMerFreqs >= repeatKmerThresholdValue);
			dynamicKmerThresholdValue = KmerThresholdTable::get(dynamicKmerSize,mode);
			fwdKmerFreqs = fwdInterval.getFreqs();
			rvcKmerFreqs = rvcInterval.getFreqs();			
			currentKmerFreqs = fwdKmerFreqs + rvcKmerFreqs;
			fixedMerFreqs = FixedMerInterval[movePos].getFreqs();
			freqsDiff = (float)fixedMerFreqs / (float)maxFixedMerFreqs;
			//Gerneral seed extension strategy. Noted by KuanWeiLee
			if  (
				   dynamicKmerSize > m_params.kmerLengthUpperBound											//1.over length
				|| isLowComplexity(kmer,GCratio)															//2.low complexity
				|| fixedMerFreqs < staticKmerThresholdValue													//3.fixed kmer frequency
				|| !OVERTHRESHOLD(currentKmerFreqs,dynamicKmerThresholdValue,fwdKmerFreqs,rvcKmerFreqs)		//4.dynamic kmer frequency
				)
			{
				if(isSeed)
				{
					kmer.erase(--dynamicKmerSize);
					seedEndPos = seedStartPos + dynamicKmerSize - 1;
					startPos = seedEndPos;
				}
				break;
			}
			//The order between general seed extentsion and kmer hitchhike strategy may make difference 
			//,which needs furher observation. Noted by KuanWeiLee 20171027
			//Kmer Hitchhike strategy.
			if(isRepeat && freqsDiff < m_params.khhRatio)													//5.hitchhiking kmer (HIGH-->LOW)
			{
				kmer.erase(--dynamicKmerSize);
				seedEndPos = seedStartPos + dynamicKmerSize - 1;
				startPos = seedEndPos + 1;
				break;
			}
			else if(fixedMerFreqs > repeatKmerThresholdValue && freqsDiff > m_params.r_khhRatio)			//6.hitchhiking kmer (LOW-->HIGH)
			{
				isSeed = false;
				startPos = movePos - 1;
				break;
			}
			
			maxFixedMerFreqs = MAX(maxFixedMerFreqs,fixedMerFreqs);
		}
		if(isSeed)
		{
			SeedFeature newSeed(seedStartPos, kmer, isRepeat, staticKmerSize, m_params.PBcoverage>>1, maxFixedMerFreqs);
			newSeed.estimateBestKmerSize(m_params.indices);
			seedVec.push_back(newSeed);
		}
	}
	seedVec = removeHitchhikingSeeds(seedVec, result);
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

PacBioSelfCorrectionProcess::SeedVector PacBioSelfCorrectionProcess::removeHitchhikingSeeds(SeedVector initSeedVec, PacBioSelfCorrectionResult& result)
{
	for(SeedVector::iterator iterQuery = initSeedVec.begin(); iterQuery != initSeedVec.end(); iterQuery++)
	{
		SeedFeature& query = *iterQuery;
		SeedVector::iterator iterTarget = iterQuery + 1;
		//if(query.isHitchhiked) continue;
		while(iterTarget != initSeedVec.end() && (*iterTarget).seedStartPos - query.seedEndPos <= m_params.repaetDistance)
		{
			SeedFeature& target = *iterTarget;
			iterTarget++;
			//if(target.isHitchhiked) continue;
			float freqsDiff = (float)target.maxFixedMerFreqs / (float)query.maxFixedMerFreqs;
			target.isHitchhiked = target.isHitchhiked || (query.isRepeat && freqsDiff < m_params.shhRatio); //HIGH --> LOW
			query.isHitchhiked = query.isHitchhiked || (target.isRepeat && freqsDiff > m_params.r_shhRatio);//LOW  --> HIGH
		}
	}
	unsigned int index = 0;
	SeedVector finalSeedVec, outcastSeedVec;
	finalSeedVec.reserve(initSeedVec.size());
	outcastSeedVec.reserve(initSeedVec.size()>>3);
	for(SeedVector::iterator iter = initSeedVec.begin(); iter != initSeedVec.end(); iter++, index++)
	{
		if((*iter).isHitchhiked)
			outcastSeedVec.push_back(*iter);
		else
			finalSeedVec.push_back(*iter);
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
	size_t seqLen = seq.length();
	size_t countG =0 ;
	size_t countC =0 ;
	size_t countT =0 ;
	size_t countA =0 ;

	for (size_t i=0; i<seqLen; i++)
	{
		switch(seq[i]){
			case 'A': countA ++ ;break;
			case 'T': countT ++ ;break;
			case 'C': countC ++ ;break;
			case 'G': countG ++ ;break;
			default:  assert(false);
		}
	}
	GCratio = (float)(countG+countC)/seqLen;
	return (float)countA/seqLen >= threshold || (float)countT/seqLen >= threshold
			|| (float)countC/seqLen >= threshold || (float)countG/seqLen >= threshold;
}

void PacBioSelfCorrectionProcess::write(std::ostream& outfile, const SeedVector& seedVec) const
{
	for(SeedVector::const_iterator iter = seedVec.begin(); iter != seedVec.end(); iter++)
		outfile << (*iter).seedStr << "\t" << (*iter).maxFixedMerFreqs << "\t" 
		<< (*iter).seedStartPos << "\t" << ((*iter).isRepeat ? "Yes" : "No") << "\n";
}
//Correct sequence by FMWalk & MSAlignment. Noted by KuanWeiLee
void PacBioSelfCorrectionProcess::initCorrect(std::string& readSeq, const SeedVector& seedVec, SeedVector& pieceVec, PacBioSelfCorrectionResult& result)
{
	std::ostream* pExtendWriter = NULL;
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
			firstFMExtensionType += 4;
			if(m_params.DebugExtend)
				*pExtendWriter << source.seedStartPos << "\t" << target.seedStartPos << "\t" << firstFMExtensionType << "\n";
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
	if(pExtendWriter != NULL) delete pExtendWriter;
}

int PacBioSelfCorrectionProcess::correctByFMExtension
(const SeedFeature& source, const SeedFeature& target, const std::string& in, std::string& out, PacBioSelfCorrectionResult& result)
{
	Timer* FMTimer = new Timer("FM Time",true);
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
	//float freqsDiff = (float)target.maxFixedMerFreqs / (float)source.maxFixedMerFreqs;
	bool isFromRtoU = false;
	//isFromRtoU = source.isRepeat && freqsDiff < m_params.shhRatio;//HIGH --> LOW
	isFromRtoU = source.isRepeat && source.maxFixedMerFreqs > target.maxFixedMerFreqs;
	if(isFromRtoU)
	{
		std::swap(src,trg);
		src = reverseComplement(src);
		trg = reverseComplement(trg);
		path = reverseComplement(path);
	}
	
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
	Timer* DPTimer = new Timer("DP Time", true);
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
	size_t totalMaxFixedMerFreqs = source.maxFixedMerFreqs + target.maxFixedMerFreqs, min_call_coverage = 15;
	identity += (totalMaxFixedMerFreqs > 50  ? 0.5 : 0);
	identity += (totalMaxFixedMerFreqs > 100 ? 0.5 : 0);
	min_call_coverage = totalMaxFixedMerFreqs > 50 ? totalMaxFixedMerFreqs * 0.4 : min_call_coverage;
	
	
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
PacBioSelfCorrectionPostProcess::PacBioSelfCorrectionPostProcess(std::string correctFile,
std::string discardFile,
const PacBioSelfCorrectionParameters params):
m_params(params),
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