#include <iostream>
#include <limits>
#include "Util.h"
#include "KmerThreshold.h"

//The threshold table is for determining correct value to tell a kmer,
//which is contributed from ChengWei Tsai & KuanWei Lee;
//there are 3 type : 0->lowcov, 1->unique, 2->repeat.
//Noted by KuanWeiLee 20180115
namespace KmerThreshold
{
	int m_startLen = 9;
	int m_endLen = 0;
	int m_coverage = 0;
	float* m_table[3] = {nullptr, nullptr, nullptr};
	float m_formula[3][6] =
	{
		//LOWCOV:0
//		{0,                0,              0,             0.05776992234, -0.4583043394, 10.19159685},//lowcov-old
//		{0,                0,              0,             0.06660714286, -0.4419117647, 10.10313375},//lowcov-new-line
		{0.0004799107143, -0.008037815126, 0.03673552754, 0.1850695903,  -1.572552521,  18.0522088 },//lowcov-new-poly
		//UNIQUE:1
//		{0,                0,              0,             0.0710704607,  -0.5445663957, 12.26253388},//unique-old
//		{0,                0,              0,             0.07928571429, -0.5568627451, 12.730007  },//unique-new-line
		{0.0002901785714, -0.009386554622, 0.04557656396, 0.252759979,   -1.997682073,  22.24817344},//unique-new-poly
//		{0.0001829519687, -0.006932900433, 0.02186460882, 0.2141542981,  -1.05755223,   13.19734651},//unique-update
		//REPEAT:2
//		{0.0007812499999, -0.1050420168,   0.218837535,   3.379382878,   -10.6797619,   139.8140931} //repeat-new-0
//		{0.004227120536,  -0.238737395,    0.7860790149,  7.077155003,   -35.45770833,  400.2657527} //repeat-new-1
//		{0.007672991071,  -0.3724327731,   1.353320495,   10.77492713,   -60.23565476,  660.7174122} //repeat-new-2
//		{0.01111886161,   -0.5061281513,   1.920561975,   14.47269925,   -85.01360119,  921.1690717} //repeat-new-3
		{0.01456473214,   -0.6398235294,   2.487803455,   18.17047138,   -109.7915476,  1181.620731} //repeat-new-4
//		{0.008677592249,  -0.4779480519,   1.437886548,   15.45189445,   -71.31698251,  839.0218299} //repeat-update
	};// x*x              x*y              y*y            x              y              (constant)
	std::ostream* pTableWriter = nullptr;
	void initialize(int start, int end, int cov, std::string dir)
	{
	//	m_startLen = start;
		m_endLen = end;
		m_coverage = cov;
		for(auto& iter : m_table)
			iter = new float[m_endLen + 1]{0};
		pTableWriter = createWriter(dir + "threshold-table");
	}
	float calculate(int type, int x, int y)
	{
		float const *formula = m_formula[type];
		return formula[0]*x*x + formula[1]*x*y + formula[2]*y*y + formula[3]*x + formula[4]*y + formula[5];
	}
	void compute()
	{
		float min_lowcov, min_unique, min_repeat;
		min_lowcov = min_unique = min_repeat = std::numeric_limits<float>::max();
		int x = m_coverage;
		for(int k = m_startLen; k <= m_endLen; k++)
		{
			int y = k;
			float lowcov = calculate(0, x, y);
			float unique = calculate(1, x, y);
			float repeat = calculate(2, x, y);
		//	lowcov = std::max(lowcov, 5);
		//	unique = std::max(unique, 5);
			lowcov = std::min(lowcov, min_lowcov);
			unique = std::min(unique, min_unique);
			repeat = std::min(repeat, min_repeat);
			m_table[0][k] = lowcov;
			m_table[1][k] = unique;
			m_table[2][k] = repeat;
			min_lowcov = std::min(lowcov, min_lowcov);
			min_unique = std::min(unique, min_unique);
			min_repeat = std::min(repeat, min_repeat);
		}
	}
	void write()
	{
		*pTableWriter << "Coverage : " << m_coverage << "\n" << "size\tlowcov\tuique\trepeat\n";
		float const *lowcov = m_table[0] + m_startLen;
		float const *unique = m_table[1] + m_startLen;
		float const *repeat = m_table[2] + m_startLen;
		for(int k = m_startLen; k <= m_endLen; k++, lowcov++, unique++, repeat++)
			*pTableWriter << k << "\t" << *lowcov << "\t" << *unique << "\t" << *repeat << "\n";
	}
	void release()
	{
		for(auto& iter : m_table)
			delete[] iter;
		delete pTableWriter;
	}
};
