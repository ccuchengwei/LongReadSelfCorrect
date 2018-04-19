#include <iostream>
#include <limits>
#include <cmath>
#include <algorithm>
#include "Util.h"
#include "KmerThreshold.h"
//The threshold table is for determining correct value to tell a kmer,
//which is contributed from ChengWei Tsai & KuanWei Lee;
//there are 3 mode : 0->lowcov, 1->unique, 2->repeat.
//Noted by KuanWeiLee 20180115
static const float formula[3][6] =
{
	//LOWCOV:0
//	{0,                0,              0,             0.05776992234, -0.4583043394, 10.19159685},//lowcov-old
	{0.0004799107143, -0.008037815126, 0.03673552754, 0.1850695903,  -1.572552521,  18.0522088 },//lowcov-new
	
	//UNIQUE:1
//	{0,                0,              0,             0.0710704607,  -0.5445663957, 12.26253388},//unique-old
//	{0.0002901785714, -0.009386554622, 0.04557656396, 0.252759979,   -1.9976820730, 22.24817344},//unique-original
	{0.0003348214286, -0.009112394958, 0.04286714686, 0.240519958,   -1.8793367350, 21.29319228},//unique-update

	//REPEAT:2
//	{0.01456473214,   -0.6398235294,   2.487803455,   18.17047138,   -109.7915476,  1181.620731} //repeat-original
	{0.01714285714,   -0.6193907563,   2.266956783,   17.28450630,   -100.6983493,  1103.571729} //repeat-update
};// x*x              x*y              y*y            x              y              (constant)

KmerThreshold::KmerThreshold():table{nullptr},pTableWriter(nullptr),set(false)
{
}

KmerThreshold::~KmerThreshold()
{
	if(pTableWriter != nullptr)
	{
		*pTableWriter << "Coverage : " << cov << "\n" << "size\tlowcov\tunique\trepeat\n";
		write(*pTableWriter);
	}
	for(auto& iter : table)
		delete iter;
	delete pTableWriter;
}

void KmerThreshold::initialize(int _start, int _end, int _cov, const std::string& _dir)
{
	if(set) return;
	set = true;
	start = std::max(_start, 15);
	end  = _end;
	cov  = _cov;
	if(!_dir.empty())
		pTableWriter = createWriter(_dir + "threshold-table");
	for(int mode = 0; mode <= 2; mode++)
	{
		table[mode] = new float[end + 2]{0};
		float* value = table[mode] + start;
		float cavity = std::numeric_limits<float>::max();
		for(int ksize = start; ksize <= end; ksize++, value++)
		{
			cavity = std::fminf(cavity, calculate(mode, cov, ksize));
			*value = cavity;
		}
	}
}

void KmerThreshold::write(std::ostream& outfile)
{
	const float* lowcov = table[0] + start;
	const float* unique = table[1] + start;
	const float* repeat = table[2] + start;
	for(int ksize = start; ksize <= end; ksize++, lowcov++, unique++, repeat++)
		outfile << ksize << "\t" << *lowcov << "\t" << *unique << "\t" << *repeat << "\n";
}

float KmerThreshold::calculate(int t, int x, int y)
{
	const float* f = formula[t];
	float v = f[0]*x*x + f[1]*x*y + f[2]*y*y + f[3]*x + f[4]*y + f[5];
	return std::fmaxf(v, 2.0f);
}
