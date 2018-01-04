// KmerFreqProcess.cpp - Check kmer distribution and error condition
//
#include "KmerFreqProcess.h"
#include "BWTAlgorithms.h"
//process
KmerFreqResult KmerFreqProcess::process(const SequenceWorkItem& workItem)
{
	KmerFreqResult result;
	std::string id = workItem.read.id, query = workItem.read.seq.toString();
	std::istream* pAlignReader = createReader(m_params.align + id + ".align");
	
	if(is_empty(pAlignReader))
	{
		delete pAlignReader;
		return result;
	}
	
	std::string key, garbage;
	int qStart, qEnd, tStart, tEnd;
	*pAlignReader
	>>	garbage	>>	key		>>	garbage
	>>	garbage	>>	garbage	>>	garbage
	>>	qStart	>>	qEnd	>>	tStart
	>>	tEnd	>>	garbage	>>	garbage;
	delete pAlignReader;
	std::string target = m_params.refMap[key];
	int	loffset = (int)((qStart - 1) * 1.3f + 0.5f) + 10,
		roffset = (int)((query.length() - qEnd) * 1.3f + 0.5f) + 10;
	bool rvc = tStart > tEnd;
	if(rvc)
	{
		std::swap(tStart,tEnd);
		std::swap(loffset,roffset);
	}
	tStart -= loffset;
	tEnd += roffset;
	tStart--;
	tEnd--;
	tStart = validatePos(tStart, target.length());
	tEnd = validatePos(tEnd, target.length());
	target = target.substr(tStart, tEnd - tStart + 1);
	if(rvc)
		target = reverseComplement(target);
	for(int i = m_params.kmerSize.first; i <= m_params.kmerSize.second; i++)
		scan(i, query, target, result);
	return result;
}
void KmerFreqProcess::scan(const int staticKmerSize, const std::string& query, const std::string& target, KmerFreqResult& result)
{
	int anchor = 0, range = (int)(staticKmerSize * 2.5f + 0.5f), section = (int)(staticKmerSize/3.f);
	for(int i = 0; i <= query.length() - staticKmerSize; i++)
	{
		std::string kmer = query.substr(i, staticKmerSize), partition;
		int kmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
		anchor = validatePos(anchor, target.length());
		partition = target.substr(anchor, range);
		int fpos = partition.find(kmer);
		if(fpos != std::string::npos)
		{
			result.correctKdMap[staticKmerSize].add(kmerFreq);
			anchor += (fpos + 1);
			range = staticKmerSize * 2;
		}
		else
		{
			if(kmerFreq == 1) continue;
			result.errorKdMap[staticKmerSize].add(kmerFreq);
			range += section;
		}
	}
}
//postprocess
KmerFreqPostProcess::KmerFreqPostProcess(KmerFreqParameters params):m_params(params)
{
	for(int i = m_params.kmerSize.first; i <= m_params.kmerSize.second; i++)
	{
		m_pCorrectWriterMap[i] = createWriter(m_params.directory + std::to_string(i) + ".correct.kf");
		m_pErrorWriterMap[i] = createWriter(m_params.directory + std::to_string(i) + ".error.kf");
		//m_pSplitCorrectWriterMap[i] = createWriter(m_params.directory + "split/" + std::to_string(i) + ".correct.kf");
		//m_pSplitErrorWriterMap[i] = createWriter(m_params.directory  + "split/" + std::to_string(i) + ".error.kf");
	}
}
KmerFreqPostProcess::~KmerFreqPostProcess()
{
	for(int i = m_params.kmerSize.first; i <= m_params.kmerSize.second; i++)
	{
		m_correctKdMap[i].write(*(m_pCorrectWriterMap[i]));
		m_errorKdMap[i].write(*(m_pErrorWriterMap[i]));
		//m_correctKdMap[i].write(*(m_pSplitCorrectWriterMap[i]), KmerDistribution::TYPE::SPLIT);
		//m_errorKdMap[i].write(*(m_pSplitErrorWriterMap[i]), KmerDistribution::TYPE::SPLIT);
	}
	for(const auto& iter : m_pCorrectWriterMap)
		delete iter.second;
	for(const auto& iter : m_pErrorWriterMap)
		delete iter.second;
	//for(const auto& iter : m_pSplitCorrectWriterMap)
	//	delete iter.second;
	//for(const auto& iter : m_pSplitErrorWriterMap)
	//	delete iter.second;
}
void KmerFreqPostProcess::process(const SequenceWorkItem& workItem, const KmerFreqResult& result)
{
	for(int i = m_params.kmerSize.first; i <= m_params.kmerSize.second; i++)
	{
		m_correctKdMap[i] += result.correctKdMap[i];
		m_errorKdMap[i] += result.errorKdMap[i];
	}
}