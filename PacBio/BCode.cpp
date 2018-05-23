#include <iostream>
#include "Util.h"
#include "BCode.h"
#define hex_num(o) (((o & 1) >> 0) + ((o & 2) >> 1) + ((o & 4) >> 2) + ((o & 8) >> 3))

const std::map<char, int> base_hex =
{
	{'a', 1}, {'t', 2}, {'c', 4}, {'g', 8},
	{'A', 1}, {'T', 2}, {'C', 4}, {'G', 8}
};

const std::map<char, int> char_int =
{
	{'0',  0}, {'1',  1}, {'2',  2}, {'3',  3},
	{'4',  4}, {'5',  5}, {'6',  6}, {'7',  7},
	{'8',  8}, {'9',  9}, {'a', 10}, {'b', 11},
	{'c', 12}, {'d', 13}, {'e', 14}, {'f', 15}
};


std::map<std::string, BCode::BCodeVector>& BCode::Log()
{
	static std::map<std::string, BCodeVector> log;
	return log;
}

void BCode::load(const std::string& barcode)
{
	if(!Log().empty())
	{
		std::cerr << "2nd loading is forbidden.\n";
		exit(EXIT_FAILURE);
	}
	std::cerr << "Loading BARCODE: " << barcode << '\n';
	std::istream* pBCodeReader = createReader(barcode);
	while(true)
	{
		if(pBCodeReader->eof()) break;
		std::string qname, tname, code, rvc, sup;
		int qstart, qend, tstart, tend;
		*pBCodeReader
		>> qname >> qstart >> qend
		>> tname >> tstart >> tend
		>> code  >> rvc    >> sup;
		Log()[qname].push_back(BCode(qstart, qend, code, (rvc == "True" ? true : false)));
	}
	delete pBCodeReader;
}

//'in[pos::step]'
std::string BCode::fetch(const std::string& in, int pos, int step)
{
	pos = getPys(pos, in.size());
	std::string out;
	for(int i = pos; (i >= 0 && i < (int)in.size()); i += step)
		out += in[i];
	return out;
}

int BCode::sum(const std::string& in)
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

//get 'py'thonic po'sition'
int BCode::getPys(int pos, int len)
{
	if(pos < 0) pos += len;
//	assert(pos >= 0 && pos < len);
	assert(pos >= 0);
	return pos;
}

bool BCode::validate(int pos, int ksize, const BCode& block, const std::string& seq)
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