#pragma once
#include "Variable.h"
#include "VariableMapping.h"

class CIndexArray
{
public:
	CIndexArray();
	~CIndexArray();
	void resetIndexArray()									{ m_numDefinedVar = 0; }
	void addVarIndex(VECTOR_ELEMENT_TYPE varIdx);
	bool isVarIndexInIndexArray(VECTOR_ELEMENT_TYPE idx) const;
private:
	void AlocateArrayForIndices();

	VECTOR_ELEMENT_TYPE *m_idxVarLast;	// indicies of the variables which was defined AND their values  
										// are axactly the same as they are defined by some equation
	size_t m_numDefinedVarMax;			// the length of array allocated for m_idxVarLast
	size_t m_numDefinedVar;				// the number of indices currently stored in m_idxVarLast
};

template<class T>
class CEquation
{
public:
	CEquation()												{ m_pIndexArray = NULL;  }
	~CEquation()											{ delete m_pIndexArray; }
	void initEquation(CVariable *pntr, T nVar, T rPart);
	void releaseEquation();
	inline CVariable *firstVariable() const					{ return m_pVar; }
	inline T numbVar() const								{ return m_nVar; }
	inline T rightPart() const								{ return m_nRightPart; }
	inline size_t adjustRightPart(T val, int mult = 1)  
															{ return setRightPart(rightPart() - mult * val); }
	int removeVariable(const T valVar, char *pVarDefined, CVariableMapping<T> *pVarValue, CEquation<T> **pDefVar, int val = 1);
	void addVariable(T valVar);
	inline bool solved() const								{ return m_bSolved; }
	size_t numDefinedVar() const							{ return m_numDefinedVar; }
	void decreaseNumDefinedVar()							{ m_numDefinedVar--; }
	bool isVarIndexInIndexArray(T i) const					{ return m_pIndexArray && m_pIndexArray->isVarIndexInIndexArray(i); }

private: 
	inline void setSolved(bool val = true)					{ m_bSolved = val; }
	inline void adjustNumbVar(T val)						{ m_nVar -= val; }
	inline void setNumbVar(T val)							{ m_nVar = val; }
	inline void setFirstVariable(CVariable *pVar)			{ m_pVar = pVar; }
	inline size_t setRightPart(T val)						{ return m_nRightPart = val; }
	inline void setNumDefinedVar(size_t value)				{ m_numDefinedVar = value; }
	void addVarIndex(T varIdx); 
	void resetIndexArray() const;

	CVariable *m_pVar;
	T m_nVar;
	T m_nRightPart;
	bool m_bSolved;
	size_t m_numDefinedVar;				// the number of variables defined from that equation
	CIndexArray *m_pIndexArray;

	MY_ID
};