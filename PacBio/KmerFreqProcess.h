// KmerFreqProcess - Check kmer distribution and error condition
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
typedef std::map<int,std::ostream*> pOstreamMap;
typedef std::map<int,KmerDistribution> kdMap;
// Parameter object
struct KmerFreqParameters
{
	KmerFreqParameters(const std::map<std::string,std::string>& seqMap):refMap(seqMap){};
	BWTIndexSet indices;
	std::string directory;
	std::string align;
	const std::map<std::string,std::string>& refMap;
	std::pair<int,int> kmerSize;
};

struct KmerFreqResult
{
	kdMap correctKdMap;
	kdMap errorKdMap;

};

//
class KmerFreqProcess
{
	public:
		KmerFreqProcess(KmerFreqParameters params):m_params(params){ };
		~KmerFreqProcess(){ };
		KmerFreqResult process(const SequenceWorkItem& workItem);
	
	private:
		void scan(int currentKmerSize, const std::string& query, const std::string& target, KmerFreqResult& result);
		inline int validatePos(int pos, int seqLen)
		{
			seqLen--;
			if(pos < 0) return 0;
			if(pos > seqLen - 1) return seqLen - 1;
			return pos;
		}
	
		KmerFreqParameters m_params;
};

//
class KmerFreqPostProcess
{
	public:
		KmerFreqPostProcess(KmerFreqParameters params);
		~KmerFreqPostProcess();

		void process(const SequenceWorkItem& workItem, const KmerFreqResult& result);

	private:
		KmerFreqParameters m_params;
		kdMap m_correctKdMap;
		kdMap m_errorKdMap;
		pOstreamMap m_pCorrectWriterMap;
		pOstreamMap m_pErrorWriterMap;
		//pOstreamMap m_pSplitCorrectWriterMap;
		//pOstreamMap m_pSplitErrorWriterMap;

};


#endif
