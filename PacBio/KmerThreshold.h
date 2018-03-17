#ifndef KMERTHRESHOLDTABLE_H
#define KMERTHRESHOLDTABLE_H

namespace KmerThreshold
{
	extern float* m_table[3];
	extern float m_formula[3][6];
	extern std::ostream* pTableWriter;
	extern int m_startLen, m_endLen, m_coverage;
	
	void initialize(int start, int end, int cov, std::string dir);
	float calculate(int type, int x, int y);
	void compute();
	void write();
	void release();
		
};

#endif