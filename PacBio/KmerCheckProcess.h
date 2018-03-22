// KmerCheckProcess - Check kmer distribution and error condition
//

#ifndef KMERFREQPROCESS_H
#define KMERFREQPROCESS_H

#include <iostream>
#include <map>
#include <utility>
#include "BWTIndexSet.h"
#include "SequenceWorkItem.h"
#include "KmerDistribution.h"
//Helpful typedef
typedef std::map<int, std::ostream*> ostreamPtrMap;
typedef std::map<int, KmerDistribution> kdMap;
// Parameter object
struct KmerCheckParameters
{
	BWTIndexSet indices;
	std::string directory;
	std::string align;
	std::map<std::string,std::string>* pSeqMap;
	std::pair<int, int> size;
};

struct KmerCheckResult
{
	kdMap correctKdMap;
	kdMap errorKdMap;
};

//
class KmerCheckProcess
{
	public:
		KmerCheckProcess(KmerCheckParameters params):m_params(params){ }
		~KmerCheckProcess(){ }
		KmerCheckResult process(const SequenceWorkItem& workItem);
	
	private:
		void scan(int currentKmerSize, const std::string& query, const std::string& target, KmerCheckResult& result);
		inline int validatePos(int pos, int seqLen)
		{
			seqLen--;
			if(pos < 0) return 0;
			if(pos > seqLen - 1) return seqLen - 1;
			return pos;
		}
	
		KmerCheckParameters m_params;
};

//
class KmerCheckPostProcess
{
	public:
		KmerCheckPostProcess(KmerCheckParameters params);
		~KmerCheckPostProcess();

		void process(const SequenceWorkItem& workItem, const KmerCheckResult& result);

	private:
		KmerCheckParameters m_params;
		kdMap m_correctKdMap;
		kdMap m_errorKdMap;
		ostreamPtrMap m_pCorrectWriterMap;
		ostreamPtrMap m_pErrorWriterMap;

};


#endif
