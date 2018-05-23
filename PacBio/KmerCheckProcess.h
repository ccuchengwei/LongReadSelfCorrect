// KmerCheckProcess - Check kmer distribution and error condition
//

#ifndef KMERCHECKPROCESS_H
#define KMERCHECKPROCESS_H

#include <iostream>
#include <map>
#include "BWTIndexSet.h"
#include "SequenceWorkItem.h"
#include "KmerDistribution.h"
#include "BCode.h"

//Helpful typedef
typedef std::map<int, KmerDistribution> kdMap;

//Parameter
struct KmerCheckParameters
{
	BWTIndexSet indices;
	std::string directory;
	int coverage;
	int lower;
	int upper;
	int step;
};

//Result
struct KmerCheckResult
{
	kdMap crtKdMap;
	kdMap errKdMap;
	std::string readid;
};

//Process
class KmerCheckProcess
{
	public:
		KmerCheckProcess(KmerCheckParameters params):m_params(params){ }
		~KmerCheckProcess(void) = default;
		KmerCheckResult process(const SequenceWorkItem& workItem);
	
	private:
		void scan(int ksize, const BCode& block, const std::string& query, KmerCheckResult& result);
		KmerCheckParameters m_params;
};

//PostProcess
class KmerCheckPostProcess
{
	public:
		KmerCheckPostProcess(KmerCheckParameters params);
		~KmerCheckPostProcess();

		void process(const SequenceWorkItem& workItem, const KmerCheckResult& result);

	private:
		KmerCheckParameters m_params;
		kdMap m_crtKdMap;
		kdMap m_errKdMap;
		std::ostream* m_pTotalWriter;
		std::ostream* m_pValueWriter;
};

#endif
