#include "Equation.h"

#if USE_EXRA_EQUATIONS
class CColOrbit;

class CColOrbInfo
{
public:
	void Init(CColOrbit *pOrb, CColOrbInfo *prev, CColOrbInfo *next)	{ m_pOrb = pOrb;  m_pPrev = prev; m_pNext = next; }
	CColOrbInfo *prev() const						{ return m_pPrev; }
	CColOrbInfo *next() const						{ return m_pNext; }
	CColOrbit *colOrb() const						{ return m_pOrb; }
private:
	CColOrbit *m_pOrb;
	CColOrbInfo *m_pPrev;
	CColOrbInfo *m_pNext;
};

typedef CContainer<CColOrbit *> CForcibleCol;

#if 0
typedef CContainer<CColOrbInfo> CColOrbInfoArray;

// Class for storage of forcibly constructed orbits of columns of the matrix
class CForcibleCol_A : public CColOrbInfoArray
{
public:
	void addForciblyConstructedColOrbit(CColOrbit *pColOrbit, S nPart, S rowNumb);
protected:
	CForcibleCol_A(size_t nVar, size_t nRow);
	~CForcibleCol_A();
	void removeForciblyConstructedColOrbit(size_t rowNumb);
private:
	inline size_t *forcibleColNumb() const					{ return m_pForcibleColNumb; }
	size_t *m_pForcibleColNumb;            // m_pForcibleColNumb[i] is the number of forcibly constructed orbits of columns for i-th row of matrix
};
#endif

typedef CArray<CEquation *, CEquation *> CEquArray;

class CEquSystem : public CForcibleCol
{
public:
	CEquSystem(size_t nVar, size_t nRow, int t);
	~CEquSystem();
	int excludeVariables(CVariableMapping *pVarValue, int adj = 1, bool exclude = true, size_t numVar = 0, const CEquArray *pEquations = NULL, bool addToAllEquation = true);
	int excludeVariables(VECTOR_ELEMENT_TYPE *pRes, VECTOR_ELEMENT_TYPE varIdx, VECTOR_ELEMENT_TYPE varValueDelta, const VECTOR_ELEMENT_TYPE *pVarMinVal, CVariableMapping *varMapping);
	int excludeVariable(VECTOR_ELEMENT_TYPE varIdx, VECTOR_ELEMENT_TYPE varValue);
	void addForcibleOrb(CColOrbit *pntr)					{ addElement(pntr); }
	void resetVariables(size_t nVar);
	inline size_t equNumb() const							{ return m_nEquNumb; }
	inline char *varDefinedPtr() const						{ return m_pVarDefined; }
	inline CVariable **varPntr() const						{ return m_ppVariable; }
	inline size_t nVar() const								{ return m_nVar; }
	bool isSolved() const;
protected:
	void allocateMemoryForEquations(size_t nEqu);
	void resetVarPtr(size_t nVar);
	void addVariable(CVariable *pVar, size_t numVar);
	void closeVarloop() const;
	void addEquation(CVariable *pFirstVar, size_t nVar, size_t lambda);
	CVariableMapping *solveExtraEquations();
private:
	inline void resetEquNumb()								{ setEquNumb(0); }
	inline void setEquNumbMax(size_t n)						{ m_nEquNumbMax = n; }
	inline void setEquNumb(size_t val)						{ m_nEquNumb = val; }
	inline size_t equNumbMax() const						{ return m_nEquNumbMax; }
	inline void setVarPntr(CVariable *p, size_t i)			{ varPntr()[i] = p; }
	inline CEquation *equation(size_t idx) const			{ return *(equPntr() + idx); }
	inline CEquation *equation(const CVariable *pVar) const	{ return *(equPntr() + *(equIndx() + pVar->equIdx())); }
	inline CEquation **equPntr() const						{ return m_pEquPntr; }
	inline VECTOR_ELEMENT_TYPE *equIndx() const				{ return m_pEquIndx; }
	inline CVariable *variable(size_t idx) const			{ return *(varPntr() + idx); }
	inline bool isVarDefined(VECTOR_ELEMENT_TYPE idx)		{ return *(varDefinedPtr() + idx) != 0; }
	inline void setVarDefined(size_t idx, char val = 1)		{ *(varDefinedPtr() + idx) = val; }
	inline CEquation *getDefEquation(size_t idx) const		{ return *(getDefEquationPntr() + idx); }
	inline void resetDefEquation(size_t idx)				{ *(getDefEquationPntr() + idx) = NULL; }
	inline CEquation **getDefEquationPntr() const			{ return m_pEquDef; }
	void releaseVariables();
	void setEquation(CEquation *pEqu, size_t idx);
	void setEquation(CEquation *pEquFrom, const CEquation *pEquTo);
	inline CEquArray *equArray() const						{ return m_pEquArray; }
	inline auto getT() const								{ return m_t; }

	CVariable **m_ppVariable;
	CEquation **m_pEquPntr;
	CEquation **m_pEquDef;         // m_pEquDef[i] is the pointer to the equation where from the value of i-th variable was defined
	CEquArray *m_pEquArray;
	VECTOR_ELEMENT_TYPE *m_pEquIndx;
	char *m_pVarDefined;
	size_t m_nVar;
	size_t m_nEquNumb;
	size_t m_nEquNumbMax;
	const uint m_t;
};

#endif
