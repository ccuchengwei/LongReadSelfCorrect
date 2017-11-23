// KmerFreqProcess.cpp - Check kmer distribution and error condition
//
#include "KmerFreqProcess.h"
#include "Timer.h"

KmerFreqResult KmerFreqProcess::process(const SequenceWorkItem& workItem)
{
	KmerFreqResult result;
	for(int i = m_params.kmerSize.first; i <= m_params.kmerSize.second; i++)
	{
		scan(i,workItem,result);
	}
	return result;
}
void KmerFreqProcess::scan(const int staticKmerSize, const SequenceWorkItem& workItem, KmerFreqResult& result)
{
	std::string readSeq = workItem.read.seq.toString();
	for(size_t i = 0; i <= readSeq.length() - staticKmerSize; i++)
	{
		std::string kmer = readSeq.substr(i,staticKmerSize);
		
	}
	/*
	result.correctKd[staticKmerSize].add(5);
	result.correctKd[staticKmerSize].add(10);
	result.errorKd[staticKmerSize].add(1);
	result.errorKd[staticKmerSize].add(2);
	*/
}
void KmerFreqPostProcess::process(const SequenceWorkItem& workItem, const KmerFreqResult& result)
{
	for(int i = m_params.kmerSize.first; i <= m_params.kmerSize.second; i++)
	{
		result.correctKd[i].write(*(m_pCorrectWriterMap[i]),KmerDistribution::TYPE::DATA);
		result.errorKd[i].write(*(m_pErrorWriterMap[i]),KmerDistribution::TYPE::DATA);
	}
}