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
//Helpful tydef
typedef std::map<int,std::ostream*> pOstreamMap;
typedef std::map<int,KmerDistribution> kdMap;
// Parameter object
struct KmerFreqParameters
{
	KmerFreqParameters(std::map<std::string,std::string>& seqMap):refMap(seqMap){};
	~KmerFreqParameters(){};
	BWTIndexSet indices;
	std::string directory;
	std::map<std::string,std::string>& refMap;
	std::pair<int,int> kmerSize;
};




struct KmerFreqResult
{
	KmerFreqResult(){};
	~KmerFreqResult(){};
	kdMap correctKd;
	kdMap errorKd;

};

//
class KmerFreqProcess
{
public:
	KmerFreqProcess(KmerFreqParameters params):m_params(params){};
	~KmerFreqProcess(){};
	KmerFreqResult process(const SequenceWorkItem& workItem);
	
private:
	void scan(int currentKmerSize, const SequenceWorkItem& workItem, KmerFreqResult& result);
	KmerFreqParameters m_params;

};

//
class KmerFreqPostProcess
{
public:
	KmerFreqPostProcess(KmerFreqParameters params, pOstreamMap& pCorrectWriterMap, pOstreamMap& pErrorWriterMap)
	:m_params(params), 
	m_pCorrectWriterMap(pCorrectWriterMap), 
	m_pErrorWriterMap(pErrorWriterMap)
	{};
	~KmerFreqPostProcess(){};

	void process(const SequenceWorkItem& workItem, const KmerFreqResult& result);

private:
	KmerFreqParameters m_params;
	pOstreamMap& m_pCorrectWriterMap;
	pOstreamMap& m_pErrorWriterMap;

};


#endif
