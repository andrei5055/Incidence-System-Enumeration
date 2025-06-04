#pragma once
#include <string> 
#include <vector>

#define DATE_TIME_TAG   "===>>"

class MatrDescr {
public:
	MatrDescr(size_t groupOrder, const char *cyclesDescr, size_t numMatrix);
	~MatrDescr()						{ delete [] cycleDescr(); }

	inline auto groupOrder() const		{ return m_groupOrder; }
	inline char* cycleDescr() const		{ return m_cyclesDescr; }
	inline auto numMatrix() const		{ return m_numMatrix; }
	inline void addNumMatrix(size_t n)	{ m_numMatrix += n; }
private:
	const size_t m_groupOrder;	// order of automorphism group
	char* m_cyclesDescr = NULL;	// cycles description
	size_t m_numMatrix;			// number of matrices with a given groupOrder/cyclesDescr
};

class MatrixDB {
public:
	~MatrixDB();
	void addMatrix(size_t groupOrder, const char *cyclesDescr, size_t numMatrix = 1);
	void transferMatrixDB(MatrixDB* matrDB);
	void addMatrixDB(const MatrixDB* matrDB);
	inline void setDescrStorage(std::vector<MatrDescr*>* pntr) {
		m_descrStorage = pntr;
	}
	inline auto descrStorage() const    { return m_descrStorage; }
protected:
	size_t reportResult(FILE* f) const;
private:
	std::vector<MatrDescr *>* m_descrStorage = NULL;
};
