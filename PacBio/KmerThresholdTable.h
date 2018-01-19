#ifndef KMERTHRESHOLDTABLE_H
#define KMERTHRESHOLDTABLE_H

/*
struct KmerThresholdTable
{
	enum TYPE{ LOWCOV, UNIQUE, REPEAT };
	typedef float* FloatPointer;
	static FloatPointer m_lowcov, m_unique, m_repeat;
	static std::ostream* pTableWriter;
	static int m_startLen, m_endLen, m_coverage;
	
	static float* get(TYPE mode);
	static void compute();
	static void write();
	static void release();
};
*/

namespace KmerThresholdTable
{
	enum TYPE{ LOWCOV = 0, UNIQUE, REPEAT };
	extern float* m_table[3];
	extern float m_formula[3][6];
	extern std::ostream* pTableWriter;
	extern int m_startLen, m_endLen, m_coverage;
	
	float calculate(int type, int x, int y);
	void compute();
	void write();
	void release();
		
};

#endif