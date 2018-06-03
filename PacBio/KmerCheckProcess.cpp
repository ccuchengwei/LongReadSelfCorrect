// KmerCheckProcess - Check kmer distribution and error condition
//
#include <algorithm>
#include <numeric>
#include "KmerCheckProcess.h"
#include "BWTAlgorithms.h"
#include "KmerFeature.h"
#include "SeedFeature.h"
#include "BCode.h"

//process
KmerCheckResult KmerCheckProcess::process(const SequenceWorkItem& workItem)
{
	KmerCheckResult result;
	std::string id  = workItem.read.id;
	std::string seq = workItem.read.seq.toString();
	result.readid = id;
	
	for(auto& iter : BCode::Log()[id])
		for(int k = m_params.lower; k <= m_params.upper; k += m_params.step)
			scan(k, iter, seq, result);
	return result;
}

void KmerCheckProcess::scan(int ksize, const BCode& block, const std::string& seq, KmerCheckResult& result)
{
	for(int pos = block.getStart(); pos <= (block.getEnd() - ksize); pos++)
	{
		KmerFeature curr(m_params.indices, seq, pos, ksize);
		assert(!curr.isFake());
		assert(curr.getFreq() != 0);
		if(curr.getFreq() == 1) continue;
		bool find = BCode::validate(pos, ksize, block, seq);
		if(find)
			result.crtKdMap[ksize].add(curr.getFreq());
		else
			result.errKdMap[ksize].add(curr.getFreq());
	}
}

//postprocess
KmerCheckPostProcess::KmerCheckPostProcess(KmerCheckParameters params):m_params(params)
{
	m_pTotalWriter = createWriter(m_params.directory + "total.box", std::ios_base::app);
	m_pValueWriter = createWriter(m_params.directory + "value.box", std::ios_base::app);
}
KmerCheckPostProcess::~KmerCheckPostProcess()
{
	for(int k = m_params.lower; k <= m_params.upper; k += m_params.step)
		compare(*m_pTotalWriter, *m_pValueWriter, m_params.coverage, k, m_crtKdMap[k], m_errKdMap[k]);
	delete m_pTotalWriter;
	delete m_pValueWriter;
}
void KmerCheckPostProcess::process(const SequenceWorkItem& workItem, const KmerCheckResult& result)
{
	for(int k = m_params.lower; k <= m_params.upper; k += m_params.step)
	{
		kdMap::const_iterator crt = result.crtKdMap.find(k);
		kdMap::const_iterator err = result.errKdMap.find(k);
		if(crt != result.crtKdMap.end()) m_crtKdMap[k] += crt->second;
		if(err != result.errKdMap.end()) m_errKdMap[k] += err->second;
	}
}
