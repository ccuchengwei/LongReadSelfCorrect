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
	
	static inline float get(int size, TYPE mode)
	{
		switch(mode)
		{
			case(LOWCOV):return m_lowcov[size];
			case(UNIQUE):return m_unique[size];
			case(REPEAT):return m_repeat[size];
			default:
				std::cout << "Wrong table type\n";
				exit(EXIT_FAILURE);
		}
	};
	
	static void compute();
	static void write();
	static void release();
};
*/

namespace KmerThresholdTable
{
	enum TYPE{ LOWCOV, UNIQUE, REPEAT };
	typedef float* FloatPointer;
	extern FloatPointer m_lowcov, m_unique, m_repeat;
	extern std::ostream* pTableWriter;
	extern int m_startLen, m_endLen, m_coverage;
	
	inline float get(int size, TYPE mode)
	{
		switch(mode)
		{
			case(LOWCOV):return m_lowcov[size];
			case(UNIQUE):return m_unique[size];
			case(REPEAT):return m_repeat[size];
			default:
				std::cout << "Wrong table type\n";
				exit(EXIT_FAILURE);
		}
	};
	
	void compute();
	void write();
	void release();
		
};

#endif