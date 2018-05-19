// CheckKmerProcess.cpp - Check kmer correctness & distribution condition
//
#include <algorithm>
#include <numeric>
#include "CheckKmerProcess.h"
#include "BWTAlgorithms.h"
#include "KmerFeature.h"
#include "SeedFeature.h"
#define hex_num(o) (((o & 1) >> 0) + ((o & 2) >> 1) + ((o & 4) >> 2) + ((o & 8) >> 3))

static const std::map<char, int> base_hex = 
{
	{'a', 1}, {'t', 2}, {'c', 4}, {'g', 8},
	{'A', 1}, {'T', 2}, {'C', 4}, {'G', 8}
};

static const std::map<char, int> char_int =
{
	{'0',  0}, {'1',  1}, {'2',  2}, {'3',  3},
	{'4',  4}, {'5',  5}, {'6',  6}, {'7',  7},
	{'8',  8}, {'9',  9}, {'a', 10}, {'b', 11},
	{'c', 12}, {'d', 13}, {'e', 14}, {'f', 15}
};
 
//process
CheckKmerResult CheckKmerProcess::process(const SequenceWorkItem& workItem)
{
	CheckKmerResult result;
	std::string id  = workItem.read.id;
	std::string seq = workItem.read.seq.toString();
	result.readid = id;
//	std::cout << id << '\n';
	const std::list<CodeBlock>& alignList = (*m_params.pAlignLog)[id];
	if(m_params.mode)
		for(auto& iter : alignList)
			for(int k = m_params.size.first; k <= m_params.size.second; k++)
				scan(k, iter, seq, result);
	else
		for(auto& e : SeedFeature::Log()[id])
		{
			int status = 2;
			for(auto& i : alignList)
			{
				if(e.seedStartPos >= i.getStart() && e.seedEndPos <= i.getEnd())
				{
					status = validate(e.seedStartPos, e.seedLen, i, seq) ? 0 : 1;
					break;
				}
			}
			result.status[status]++;
		}
	return result;
}

void CheckKmerProcess::scan(int ksize, const CodeBlock& block, const std::string& seq, CheckKmerResult& result)
{
	for(int pos = block.getStart(); pos <= (block.getEnd() - ksize); pos++)
	{
		KmerFeature curr(m_params.indices, seq, pos, ksize);
		assert(!curr.getPseudo());
		assert(curr.getFreq() != 0);
		if(curr.getFreq() == 1) continue;
		bool find = validate(pos, ksize, block, seq);
		if(find)
			result.corKdMap[ksize].add(curr.getFreq());
		else
			result.errKdMap[ksize].add(curr.getFreq());
//		std::cout << pos << '\t' << curr.getWord() << '\t' << (find ? "True\n" : "False\n");
	}
}

std::string CheckKmerProcess::fetch(const std::string& in, int pos, int step)
{
	pos = getPys(pos, in.size());
	std::string out;
	for(int i = pos; (i >= 0 && i < (int)in.size()); i += step)
		out += in[i];
	return out;
}

int CheckKmerProcess::sum(const std::string& in)
{
	int out = 0;
	for(char c : in)
	{
	//	int v = char_int[c];
		int v = char_int.at(c);
		out += v;
	}
	return out;
}

int CheckKmerProcess::getPys(int pos, int len)
{
	if(pos < 0) pos += len;
//	assert(pos >= 0 && pos < len);
	assert(pos >= 0);
	return pos;
}

bool CheckKmerProcess::validate(int pos, int ksize, const CodeBlock& block, const std::string& seq)
{
	int start = pos;
	int end   = start + ksize;
	int base  = block.getStart();
	int first = (start - base) * 2;
	int last  = (end   - base) * 2 - 1;
	const std::string kmer = seq.substr(pos, ksize);
	const std::string code = block.getCode();
	const std::string info = code.substr(first, (last - first));
//	std::cout << first + 1 << '-' << last << '\n' << info << '\n';
	bool rvc = block.getRvc();
	int sign = rvc ?   - 1 :   1; // sign -> fwd :   1; rvc :    -1
	int bit  = rvc ?     0 :   1; // bit  -> fwd :   1; rvc :     0
	int pole = rvc ? start : end; // pole -> fwd : end; rvc : start
//	INSERTION GAP
	int upper = sum(fetch(info, 0, 2));
	if(upper > 0)
	{
		int igap  = 0;
		int n     = 0;
		for(char c : fetch(info, -bit, -sign*2))
		{
		//	int v = char_int[c];
			int v = char_int.at(c);
			if( ! ( (igap == 0 && (v == 0 || v == 1)) || (igap > 0 && v == 1) ) ) break;
			n += 1;
			igap += v;
		}
		if((upper - igap) != 0) return false;
		if(igap > 0)
		{
			int ioffset = 0;
			for(char c : fetch(fetch(code, 0, 2), (pole - base + bit - 1), sign))
			{
			//	int v = char_int[c];
				int v = char_int.at(c);
				if(v != 1) break;
				ioffset += 1;
			}
			if((n - igap) > 0 && ioffset > 0) return false;
			for(int i = 0; i < n; i++)
			{
				if	(!(
						fetch(code, 0, 2)[pole - base + sign*(1 - bit + ioffset + i) - sign*(n - igap)] == '0'
					&& 	kmer[getPys(-sign*(n + bit - 1 - i), ksize)] == seq[pole + sign*(1 - bit + ioffset + i) - sign*(n - igap)]
					) ) return false;
			}
		}
	}
//	DELETION GAP
	int lower = sum(fetch(info, 1, 2));
	if(lower > 0)
	{
		int dgap  = 0;
		int m     = 0;
		int hex   = 0;
		for(char c : fetch(info, -sign*(1 + bit), -sign*2))
		{
		//	int v = char_int[c];
			int v = char_int.at(c);
			if(dgap != 0) break;
		//	hex = hex | base_hex[kmer[getPys(-sign*(bit + m), ksize)]];
			hex = hex | base_hex.at(kmer[getPys(-sign*(bit + m), ksize)]);
			m += 1;
			dgap += v;
		}
		if((lower - dgap) != 0) return false;
		if(dgap > 0)
			if( ! ( dgap == hex || (m == 1 && (dgap & hex) > 0 && hex_num(dgap) == 2 ) ) ) return false;
	}
//	PASS ALL CHECK
	return true;
}

//postprocess
CheckKmerPostProcess::CheckKmerPostProcess(CheckKmerParameters params):m_params(params)
{
	if(m_params.mode)
		for(int k = m_params.size.first; k <= m_params.size.second; k++)
		{
			m_pCorWriterMap[k] = createWriter(m_params.directory + std::to_string(k) + ".cor.kf");
			m_pErrWriterMap[k] = createWriter(m_params.directory + std::to_string(k) + ".err.kf");
		}
	else
	{
		m_status = new int[3]{0};
		m_pStatusWriter = fopen((m_params.directory + "total.stat").c_str(), "w");
	}
}
CheckKmerPostProcess::~CheckKmerPostProcess()
{
	if(m_params.mode)
	{
		for(int k = m_params.size.first; k <= m_params.size.second; k++)
		{
			m_corKdMap[k].write(*(m_pCorWriterMap[k]));
			m_errKdMap[k].write(*(m_pErrWriterMap[k]));
		}
		for(const auto& iter : m_pCorWriterMap)
			delete iter.second;
		for(const auto& iter : m_pErrWriterMap)
			delete iter.second;
	}
	else
	{
		int sum = std::accumulate(m_status, (m_status + 3), 0);
		float cor = (float)(100*m_status[0])/sum;
		float err = (float)(100*m_status[1])/sum;
		float non = (float)(100*m_status[2])/sum;
		printf("TOTAL(%d) %.2f%% %.2f%% %.2f%%\n", sum, cor, err, non);
		delete m_status;
		fclose(m_pStatusWriter);
	}
}
void CheckKmerPostProcess::process(const SequenceWorkItem& workItem, const CheckKmerResult& result)
{
	if(m_params.mode)
		for(int k = m_params.size.first; k <= m_params.size.second; k++)
		{
			kdMap::const_iterator cor = result.corKdMap.find(k);
			kdMap::const_iterator err = result.errKdMap.find(k);
			if(cor != result.corKdMap.end()) m_corKdMap[k] += cor->second;
			if(err != result.errKdMap.end()) m_errKdMap[k] += err->second;
		//	m_corKdMap[k] += result.corKdMap[k];
		//	m_errKdMap[k] += result.errKdMap[k];
		}
	else
	{
		int sum = std::accumulate(result.status, (result.status + 3), 0);
		float cor = (float)(100*result.status[0])/sum;
		float err = (float)(100*result.status[1])/sum;
		float non = (float)(100*result.status[2])/sum;
		if(result.status[1] > 0)
			fprintf(m_pStatusWriter, "%s %.2f%% %.2f%% %.2f%%\n", result.readid.c_str(), cor, err, non);
	//	for(int i = 0; i < 3; i++)
	//		m_status[i] += result.status[i];
		std::transform(m_status, (m_status + 3), result.status, m_status, std::plus<int>());
	}
}
