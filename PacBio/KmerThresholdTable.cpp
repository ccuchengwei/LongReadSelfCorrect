#include <iostream>
#include <limits>
#include "KmerThresholdTable.h"
#define MIN(a,b) ( (a < b) ? a : b)
/*
int KmerThresholdTable::m_startLen = 0;
int KmerThresholdTable::m_endLen = 0;
int KmerThresholdTable::m_coverage = 0;
float* KmerThresholdTable::m_lowcov = NULL;
float* KmerThresholdTable::m_unique = NULL;
float* KmerThresholdTable::m_repeat = NULL;
std::ostream* KmerThresholdTable::pTableWriter = NULL;

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

namespace KmerThresholdTable
{
	int m_startLen = 0;
	int m_endLen = 0;
	int m_coverage = 0;
	FloatPointer m_lowcov = NULL;
	FloatPointer m_unique = NULL;
	FloatPointer m_repeat = NULL;
	std::ostream* pTableWriter = NULL;
	
	float* get(TYPE mode)
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
	void compute()
	{
		float min_unique, min_repeat;
		min_unique = min_repeat = std::numeric_limits<float>::max();
		int x = m_coverage;
		for(int k = m_startLen; k <= m_endLen; k++)
		{
			int y = k;
			float lowcov = 0.05776992234f*x - 0.4583043394f*y + 10.19159685f;
			lowcov = lowcov < 5 ? 5 : lowcov;
			float unique = 0.000234375f*x*x - 0.009113445378f*x*y + 0.04496381886*y*y + 0.2529766282*x - 1.98467437*y + 22.10684816;
			//float unique = 0.0710704607f*x - 0.5445663957f*y + 12.26253388f;
			//unique = unique < 5 ? 5 : unique;
			float repeat = 0.007533482143f*x*x - 0.2664117647f*x*y + 1.200805322f*y*y + 7.283483456f*x - 59.01653361*y + 763.5592525;
			unique = MIN(unique,min_unique);
			repeat = MIN(repeat,min_repeat);
			m_lowcov[k] = lowcov;
			m_unique[k] = unique;
			m_repeat[k] = repeat;
			min_unique = MIN(unique,min_unique);
			min_repeat = MIN(repeat,min_repeat);
		}
	}
	void write()
	{
		//for(int i = m_startLen; i <= m_endLen; i++)
		//	std::cout << i << "\t" << m_lowcov[i] << "\t" << m_unique[i] << "\t" << m_repeat[i] << "\n";
		
		*pTableWriter << "Coverage : " << m_coverage << "\n" << "size\tlowcov\tuique\trepeat\n";
		float const *lowcov = m_lowcov + m_startLen;
		float const *unique = m_unique + m_startLen;
		float const *repeat = m_repeat + m_startLen;
		for(int k = m_startLen; k <= m_endLen; k++, lowcov++, unique++, repeat++)
			*pTableWriter << k << "\t" << *lowcov << "\t" << *unique << "\t" << *repeat << "\n";
		
	}
	void release()
	{
		delete[] m_lowcov;
		delete[] m_unique;
		delete[] m_repeat;
		delete pTableWriter;
	}
};
