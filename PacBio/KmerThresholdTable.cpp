#include <iostream>
#include <limits>
#include "KmerThresholdTable.h"
#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a > b ? a : b)
/*
int KmerThresholdTable::m_startLen = 0;
int KmerThresholdTable::m_endLen = 0;
int KmerThresholdTable::m_coverage = 0;
float* KmerThresholdTable::m_lowcov = nullptr;
float* KmerThresholdTable::m_unique = nullptr;
float* KmerThresholdTable::m_repeat = nullptr;
std::ostream* KmerThresholdTable::pTableWriter = nullptr;

float* KmerThresholdTable::get(TYPE mode)
{
	switch(mode)
	{
		case LOWCOV: return m_lowcov;
		case UNIQUE: return m_unique;
		case REPEAT: return m_repeat;
		default:
			std::cout << "Wrong table type\n";
			exit(EXIT_FAILURE);
	}
}
void KmerThresholdTable::compute()
{
	int x = m_coverage;
	for(int k = m_startLen; k <= m_endLen; k++)
	{
		int y = k;
		float lowcov = 0.05776992234f*x - 0.4583043394f*y + 10.19159685f;
		lowcov = lowcov < 5 ? 5 : lowcov;
		float unique = 0.000234375f*x*x - 0.009113445378f*x*y + 0.04496381886*y*y + 0.2529766282*x - 1.98467437*y + 22.10684816;
		float repeat = 0.007533482143f*x*x - 0.2664117647f*x*y + 1.200805322f*y*y + 7.283483456f*x - 59.01653361*y + 763.5592525;
		m_lowcov[k] = lowcov;
		m_unique[k] = unique;
		m_repeat[k] = repeat;
	}
}
void KmerThresholdTable::write()
{
	*pTableWriter << "Coverage : " << m_coverage << "\n" ;
	//FloatPointer lowcov , unique ,repeat;
	float const *lowcov = m_lowcov + m_startLen;
	float const *unique = m_unique + m_startLen;
	float const *repeat = m_repeat + m_startLen;
	for(int k = m_startLen; k <= m_endLen; k++, lowcov++, unique++, repeat++)
		*pTableWriter << k << "\t" << *lowcov << "\t" << *unique << "\t" << *repeat << "\n";
}
void KmerThresholdTable::release()
{
	delete[] m_lowcov;
	delete[] m_unique;
	delete[] m_repeat;
	delete pTableWriter;
}
*/
//The threshold table is for determining correct value to tell a kmer,
//which is contributed from ChengWei Tsai & KuanWei Lee
//there are 3 type : 0->lowcov, 1->unique, 2->repeat.
//Noted by KuanWeiLee 20180115
namespace KmerThresholdTable
{
	int m_startLen = 9;
	int m_endLen = 0;
	int m_coverage = 0;
	float* m_table[3] = {nullptr, nullptr, nullptr};
	float m_formula[3][6] =
	{
		//LOWCOV:0
		{0,                0,              0,             0.05776992234, -0.4583043394, 10.19159685},//lowcov-old
		//UNIQUE:1
//		{0,                0,              0,             0.0710704607,  -0.5445663957, 12.26253388},//unique-old
//		{0,                0,              0,             0.07928571429, -0.5568627451, 12.730007  },//unique-new-line
		{0.0002901785714, -0.009386554622, 0.04557656396, 0.252759979,   -1.997682073,  22.24817344},//unique-new-poly
		//REPEAT:2
//		{0.0007812499999, -0.1050420168,   0.218837535,   3.379382878,   -10.6797619,   139.8140931} //repeat-new-0
//		{0.004227120536,  -0.238737395,    0.7860790149,  7.077155003,   -35.45770833,  400.2657527} //repeat-new-1
//		{0.007672991071,  -0.3724327731,   1.353320495,   10.77492713,   -60.23565476,  660.7174122} //repeat-new-2
//		{0.01111886161,   -0.5061281513,   1.920561975,   14.47269925,   -85.01360119,  921.1690717} //repeat-new-3
		{0.01456473214,   -0.6398235294,   2.487803455,   18.17047138,   -109.7915476,  1181.620731} //repeat-new-4
	};// x*x              x*y              y*y            x              y              (constant)
	std::ostream* pTableWriter = nullptr;
	
	float calculate(int type, int x, int y)
	{
		float const *formula = m_formula[type];
		return formula[0]*x*x + formula[1]*x*y + formula[2]*y*y + formula[3]*x + formula[4]*y + formula[5];
	}
	void compute()
	{
		float min_unique, min_repeat;
		min_unique = min_repeat = std::numeric_limits<float>::max();
		int x = m_coverage;
		for(int k = m_startLen; k <= m_endLen; k++)
		{
			int y = k;
			float lowcov = calculate(0, x, y);
			float unique = calculate(1, x, y);
			float repeat = calculate(2, x, y);
			lowcov = MAX(lowcov, 5);
			//unique = MAX(unique, 5);
			unique = MIN(unique, min_unique);
			repeat = MIN(repeat, min_repeat);
			m_table[0][k] = lowcov;
			m_table[1][k] = unique;
			m_table[2][k] = repeat;
			min_unique = MIN(unique, min_unique);
			min_repeat = MIN(repeat, min_repeat);
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
