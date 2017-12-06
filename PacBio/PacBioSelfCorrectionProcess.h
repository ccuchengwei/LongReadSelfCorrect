//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioSelfCorrectionProcess - Self-correction using FM-index walk for PacBio reads
//

#ifndef PACBIOSELFCORRECTIONPROCESS_H
#define PACBIOSELFCORRECTIONPROCESS_H

#include "HashMap.h"
#include "Util.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "Metrics.h"
#include "BWTIndexSet.h"
#include "SampledSuffixArray.h"
#include "BWTAlgorithms.h"
#include "KmerDistribution.h"
#include "LongReadCorrectByOverlap.h"
#include "SeedFeature.h"

// Parameter object for the error corrector
struct PacBioSelfCorrectionParameters
{
	//PacBioSelfCorrectionAlgorithm algorithm;
	BWTIndexSet indices;

	int numKmerRounds;
	unsigned int kmerLength;
	unsigned int kmerLengthUpperBound;
	unsigned int repaetDistance = 100;//70-->100. Need further check. Noted by KuanWeiLee 11/25
	float hhRatio = 0.6f;
	float r_hhRatio = 1.f/hhRatio;
	//KmerThreshold table;

	// tree search parameters
	int maxLeaves;
	int minOverlap;
	int maxOverlap;
    int idmerLength;
    

	// PACBIO
	int minKmerLength;
	int FMWKmerThreshold;
	int seedKmerThreshold;
	int numOfNextTarget;
	int collectedSeeds;
    double ErrorRate;
	bool isSplit;
	bool isFirst;
	size_t maxSeedInterval;
	int PBcoverage;
    
    bool DebugExtend;
    bool DebugSeed;
	bool OnlySeed;
	bool NoDp;
	std::string directory;
};




struct PacBioSelfCorrectionResult
{
	PacBioSelfCorrectionResult()
	: merge(false),
	totalReadsLen(0),
	correctedLen(0),
	totalSeedNum(0),
	totalWalkNum(0),
	correctedNum(0),
	highErrorNum(0),
	exceedDepthNum(0),
	exceedLeaveNum(0),
    DPNum(0),
	seedDis(0),
    Timer_Seed(0),
    Timer_FM(0),
    Timer_DP(0) {}

	std::string readid;
	KmerDistribution kd;
	
	bool merge;
	
	size_t kmerLength;

	// PacBio reads correction by Ya, v20151001.
	std::vector<DNAString> correctedPacbioStrs;
	int64_t totalReadsLen;
	int64_t correctedLen;
	int64_t totalSeedNum;
	int64_t totalWalkNum;
	int64_t correctedNum;
	int64_t highErrorNum;
	int64_t exceedDepthNum;
	int64_t exceedLeaveNum;
    int64_t DPNum;
	int64_t seedDis;
    double Timer_Seed;
    double Timer_FM;
    double Timer_DP;
};

//
class PacBioSelfCorrectionProcess
{
public:
	PacBioSelfCorrectionProcess(const PacBioSelfCorrectionParameters params):m_params(params){};
	~PacBioSelfCorrectionProcess(){};
	PacBioSelfCorrectionResult process(const SequenceWorkItem& workItem);

private:
    FMextendParameters FMextendParameter();
    

    std::vector<SeedFeature> hybridSeedingFromPB(const std::string& readSeq, PacBioSelfCorrectionResult &result);
	void initCorrect(std::string& readSeq, std::vector<SeedFeature>& seeds, std::vector<SeedFeature>& pacbioCorrectedStrs, PacBioSelfCorrectionResult& result);
	// Perform FMindex extension between source and target seeds
	// Return FMWalkReturnType
	int extendBetweenSeeds(SeedFeature& source, SeedFeature& target, std::string& rawSeq, std::string& mergedseq,
							size_t smallKmerSize, size_t dis_between_src_target,FMextendParameters FMextendParameter, PacBioSelfCorrectionResult& result);	
	
	//bool isCloseTorepeat(std::vector<BWTIntervalPair> FixedMerInterval,size_t &currpos);
	// kmers around repeat seeds are often error seeds, split the repeat regions into high-confident seeds
	// return kmer freq of beginning and ending kmers
	//std::pair<size_t, size_t> refineRepeatSeed(const std::string readSeq, size_t& seedStartPos, size_t& seedEndPos,size_t normal_freqs);
	// return complexity of seq, default: 0.9
	bool  isLowComplexity (const std::string& seq, float & GCratio, float threshold=0.7);

	// return <0: give up and break
	// return 0: retry the same target
	// return >0: continue to next target
	//int  FMWalkFailedActions (size_t& smallKmerSize, size_t& numOfTrials, SeedFeature& sourceStrLen, SeedFeature& target, int FMWalkReturnType, int prevFMWalkReturnType);

                            
	PacBioSelfCorrectionParameters m_params;
	//std::pair<size_t,size_t> alnscore;
	
};

//postprocess
class PacBioSelfCorrectionPostProcess
{
public:
	PacBioSelfCorrectionPostProcess(std::string correctFile, std::string discardFile, const PacBioSelfCorrectionParameters params);
	~PacBioSelfCorrectionPostProcess();
	void process(const SequenceWorkItem& item, const PacBioSelfCorrectionResult& result);
	
private:

	std::ostream* m_pCorrectWriter;
	std::ostream* m_pDiscardWriter;
	std::ostream* m_pKdWriter;
	PacBioSelfCorrectionParameters m_params;
	KmerDistribution m_kd;
	
	int64_t m_totalReadsLen;
	int64_t m_correctedLen;
	int64_t m_totalSeedNum;
	int64_t m_totalWalkNum;
	int64_t m_correctedNum;
	int64_t m_highErrorNum;
	int64_t m_exceedDepthNum;
	int64_t m_exceedLeaveNum;
    int64_t m_DPNum;
	int64_t m_seedDis;
    double m_Timer_Seed;
    double m_Timer_FM;
    double m_Timer_DP;

};


#endif
