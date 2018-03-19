#include <iostream>
#include <limits>
#include "Util.h"
#include "KmerThreshold.h"
//The threshold table is for determining correct value to tell a kmer,
//which is contributed from ChengWei Tsai & KuanWei Lee;
//there are 3 type : 0->lowcov, 1->unique, 2->repeat.
//Noted by KuanWeiLee 20180115
static const float formula[3][6] =
{
	//LOWCOV:0
//	{0,                0,              0,             0.05776992234, -0.4583043394, 10.19159685},//lowcov-old
//	{0,                0,              0,             0.06660714286, -0.4419117647, 10.10313375},//lowcov-new-line
	{0.0004799107143, -0.008037815126, 0.03673552754, 0.1850695903,  -1.572552521,  18.0522088 },//lowcov-new-poly
	
	//UNIQUE:1
//	{0,                0,              0,             0.0710704607,  -0.5445663957, 12.26253388},//unique-old
//	{0,                0,              0,             0.07928571429, -0.5568627451, 12.730007  },//unique-new-line
	{0.0002901785714, -0.009386554622, 0.04557656396, 0.252759979,   -1.997682073,  22.24817344},//unique-new-poly
//	{0.0001829519687, -0.006932900433, 0.02186460882, 0.2141542981,  -1.05755223,   13.19734651},//unique-update

	//REPEAT:2
//	{0.0007812499999, -0.1050420168,   0.218837535,   3.379382878,   -10.6797619,   139.8140931} //repeat-new-0
//	{0.004227120536,  -0.238737395,    0.7860790149,  7.077155003,   -35.45770833,  400.2657527} //repeat-new-1
//	{0.007672991071,  -0.3724327731,   1.353320495,   10.77492713,   -60.23565476,  660.7174122} //repeat-new-2
//	{0.01111886161,   -0.5061281513,   1.920561975,   14.47269925,   -85.01360119,  921.1690717} //repeat-new-3
	{0.01456473214,   -0.6398235294,   2.487803455,   18.17047138,   -109.7915476,  1181.620731} //repeat-new-4
//	{0.008677592249,  -0.4779480519,   1.437886548,   15.45189445,   -71.31698251,  839.0218299} //repeat-update
};// x*x              x*y              y*y            x              y              (constant)

KmerThreshold::KmerThreshold():table{nullptr},pTableWriter(nullptr)
{
}

KmerThreshold::~KmerThreshold()
{
	*pTableWriter << "Coverage : " << cov << "\n" << "size\tlowcov\tuique\trepeat\n";
	float const *lowcov = table[0] + start;
	float const *unique = table[1] + start;
	float const *repeat = table[2] + start;
	for(int k = start; k <= end; k++, lowcov++, unique++, repeat++)
		*pTableWriter << k << "\t" << *lowcov << "\t" << *unique << "\t" << *repeat << "\n";
	
	for(auto& iter : table)
		delete iter;
	delete pTableWriter;
}

void KmerThreshold::set(int _start, int _end, int _cov, std::string& dir)
{
	if(pTableWriter != nullptr) return;
//	start = _start;
	start = 15;
	end  = _end;
	cov  = _cov;
	pTableWriter = createWriter(dir + "threshold-table");
	for(int type = 0; type <= 2; type++)
	{
		table[type] = new float[end + 2]{0};
		float cavity = std::numeric_limits<float>::max();
		for(int ksize = start; ksize <= end; ksize++)
		{
			cavity = std::min(cavity, calculate(type, cov, ksize));
			table[type][ksize] = cavity;
		}
	}
}

void KmerThreshold::print()
{
	float const *lowcov = table[0] + start;
	float const *unique = table[1] + start;
	float const *repeat = table[2] + start;
	for(int k = start; k <= end; k++, lowcov++, unique++, repeat++)
		std::cout << k << "\t" << *lowcov << "\t" << *unique << "\t" << *repeat << "\n";
}

float KmerThreshold::calculate(int t, int x, int y)
{
	float const *f = formula[t];
	return f[0]*x*x + f[1]*x*y + f[2]*y*y + f[3]*x + f[4]*y + f[5];
}
