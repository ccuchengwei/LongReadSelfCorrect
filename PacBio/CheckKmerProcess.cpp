// CheckKmerProcess.cpp - Check kmer correctness & distribution condition
//
#include "CheckKmerProcess.h"
#include "BWTAlgorithms.h"
#include "KmerFeature.h"
#include <algorithm>
#define hex_num(o) (((o & 1) >> 0) + ((o & 2) >> 1) + ((o & 4) >> 2) + ((o & 8) >> 3))

const static std::map<char, int> base_hex = 
{
	{'a', 1}, {'t', 2}, {'c', 4}, {'g', 8},
	{'A', 1}, {'T', 2}, {'C', 4}, {'G', 8}
};

const static std::map<char, int> char_int =
{
	{'0', 0},  {'1', 1},  {'2', 2},  {'3', 3},
	{'4', 4},  {'5', 5},  {'6', 6},  {'7', 7},
	{'8', 8},  {'9', 9},  {'a', 10}, {'b', 11},
	{'c', 12}, {'d', 13}, {'e', 14}, {'f', 15}
};
 
//process
CheckKmerResult CheckKmerProcess::process(const SequenceWorkItem& workItem)
{
	CheckKmerResult result;
	std::string id    = workItem.read.id;
	std::string seq = workItem.read.seq.toString();
//	std::cout << id << '\n';
	const std::list<CodeBlock>& alignList = (*m_params.pAlignRec)[id];
	for(auto& iter : alignList)
		for(int k = m_params.size.first; k <= m_params.size.second; k++)
			scan(k, iter, seq, result);
	return result;
}

void CheckKmerProcess::scan(int ksize, const CodeBlock& block, const std::string& seq, CheckKmerResult& result)
{
	for(int pos = block.getStart(); pos <= (block.getEnd() - ksize); pos++)
	{
		KmerFeature curr(m_params.indices, seq, pos, ksize);
		assert(!curr.getPseudo());
		if(curr.getFreq() == 1) continue;
		bool find = validate(pos, ksize, block, seq);
		if(find)
			result.correctKdMap[ksize].add(curr.getFreq());
		else
			result.errorKdMap[ksize].add(curr.getFreq());
//		std::cout << pos << '\t' << curr.getWord() << '\t' << (find ? "True\n" : "False\n");
	}
}

std::string CheckKmerProcess::fetch(const std::string& in, int pos, int step)
{
	pos = getPys(pos, in.size());
	std::string out;
	for(int i = pos; (i >= 0 && i < in.size()); i += step)
		out += in[i];
	return out;
}

int CheckKmerProcess::sum(const std::string& in)
{
	int out = 0;
	for(char c : in)
	{
		int v = char_int[c];
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
			int v = char_int[c];
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
				int v = char_int[c];
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
			int v = char_int[c];
			if(dgap != 0) break;
			hex = hex | base_hex[kmer[getPys(-sign*(bit + m), ksize)]];
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
	for(int k = m_params.size.first; k <= m_params.size.second; k++)
	{
		m_pCorrectWriterMap[k] = createWriter(m_params.directory + std::to_string(k) + ".correct.kf");
		m_pErrorWriterMap[k]   = createWriter(m_params.directory + std::to_string(k) + ".error.kf");
	}
}
CheckKmerPostProcess::~CheckKmerPostProcess()
{
	for(int k = m_params.size.first; k <= m_params.size.second; k++)
	{
		m_correctKdMap[k].write(*(m_pCorrectWriterMap[k]));
		m_errorKdMap[k]  .write(*(m_pErrorWriterMap[k]));
	}
	for(const auto& iter : m_pCorrectWriterMap)
		delete iter.second;
	for(const auto& iter : m_pErrorWriterMap)
		delete iter.second;
}
void CheckKmerPostProcess::process(const SequenceWorkItem& workItem, const CheckKmerResult& result)
{
	for(int k = m_params.size.first; k <= m_params.size.second; k++)
	{
		m_correctKdMap[k] += result.correctKdMap[k];
		m_errorKdMap[k]   += result.errorKdMap[k];
	}
}
