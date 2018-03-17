
#ifndef LONGREADPROBE_H
#define LONGREADPROBE_H

#include "SeedFeature.h"

struct ProbeParameters
{
	ProbeParameters(
			BWTIndexSet _indices,
			std::string _readid,
			std::string _directory,
			int _startKmerLen,
			int _scanKmerLen,
			int _kmerLenUpBound,
			int _PBcoverage,
			int _mode,
			int _repeatDis,
			float _hhRatio,
			std::array<int, 3> _kmerOffset,
			std::set<int> _kmerPool,
			bool _DebugSeed,
			bool _Manual)
	:	indices(_indices),
		readid(_readid),
		directory(_directory),
		startKmerLen(_startKmerLen),
		scanKmerLen(_scanKmerLen),
		kmerLenUpBound(_kmerLenUpBound),
		PBcoverage(_PBcoverage),
		mode(_mode),
		repeatDis(_repeatDis),
		hhRatio(_hhRatio),
		kmerOffset(_kmerOffset),
		kmerPool(_kmerPool),
		DebugSeed(_DebugSeed),
		Manual(_Manual){ }
		
	~ProbeParameters(void) = default;
		
	BWTIndexSet indices;
	std::string readid;
	std::string directory;
	int startKmerLen;
	int scanKmerLen;
	int kmerLenUpBound;
	int PBcoverage;
	int mode;
	int repeatDis;
	float hhRatio;
	std::array<int, 3> kmerOffset;
	std::set<int> kmerPool;
	bool DebugSeed;
	bool Manual;
};

namespace LongReadProbe
{	
	void searchSeedsWithHybridKmers(const std::string& readSeq, SeedFeature::SeedVector& seedVec, const ProbeParameters& m_params);
	void getSeqAttribute(const std::string& seq, int* const type, const ProbeParameters& m_params);
	SeedFeature::SeedVector removeHitchhikingSeeds(SeedFeature::SeedVector initSeedVec, int const *type, const ProbeParameters& m_params);
	void write(std::ostream& outfile, const SeedFeature::SeedVector& seedVec);
};
#endif
