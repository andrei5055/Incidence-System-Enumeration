#include "matrix.h"
#include "RowSolution.h"
#include "InSysSolver.h"
#include "ColOrbitManager.h"
#include "CanonicityChecker.h"
#include "ThreadEnumerator.h"

class CEnumerator : public CColOrbitManager
{
public:
	CEnumerator(const CMatrix *pMatrix, bool matrOwner = false);
	virtual ~CEnumerator();
	inline CRowSolution *rowStuff(size_t nRow = 0) const    { return m_pRow[nRow]; }
	ulonglong Enumerate(bool writeFile = false, CEnumInfo *pEnumInfo = NULL, const CEnumerator *pMaster = NULL, t_threadCode *pTreadCode = NULL);
	inline const CMatrix *matrix() const                    { return m_pMatrix; }
    virtual bool makeJobTitle(char *buffer, int len, const char *comment = "") const	
															{ return false; }
	virtual VECTOR_ELEMENT_TYPE getX0_3() const             { return 0; }
	inline CRowSolution **rowStuffPntr() const				{ return m_pRow;  }
	inline void closeFile()									{ if (outFile()) { fclose(outFile()); setOutFile(NULL); }}
    virtual size_t firstUnforcedRow() const                 { return 0; }
	virtual void setFirstUnforcedRow(size_t rowNum = 0)     {}
	virtual size_t *forcibleLambdaPntr() const              { return NULL; }
	virtual bool isIS_enumerator() const                    { return false; }
	inline void setMatrOwner(bool val)						{ m_bMatrOwner = val; }
    inline CCanonicityChecker *canonChecker() const         { return m_pCanonChecker; }
protected:
	inline C_InSys *getInSys() const                        { return isIS_enumerator()? (C_InSys *)matrix() : NULL; }
    virtual void setX0_3(VECTOR_ELEMENT_TYPE value)         {}
	inline CRowEquation *rowEquation() const                { return m_pRowEquation; }
	inline void initRowEquation(unsigned int len)           { setRowEquation(new CRowEquation(len)); }
	virtual int unforcedElement(const CColOrbit *p, int nRow) const    { return -1; }
    virtual bool sortSolutions(CRowSolution *p, size_t idx) { return false;  /* not implemented */ }
	inline void setRowEquation(CRowEquation *pntr)          { m_pRowEquation = pntr; }
    inline size_t rowNumb() const                           { return matrix()->rowNumb(); }
    inline void setCanonChecker(CCanonicityChecker *p)      { m_pCanonChecker = p; }
	virtual bool makeFileName(char *buffer, int len, const char *ext = NULL) const	{ return false; }
	inline void setUseCanonGroup(bool val)					{ m_bUseCanogGroup = val; }
	inline bool useCanonGroup() const						{ return m_bUseCanogGroup; }
	virtual bool noReplicatedBlocks() const					{ return false; }
	CColOrbit *MakeRow(const VECTOR_ELEMENT_TYPE *pRowSolution) const;
private:
	virtual CRowSolution *setFirstRowSolutions()            { return NULL; }
	CRowSolution *FindRowSolution(PERMUT_ELEMENT_TYPE lastRightPartIndex = -1);
	virtual size_t MakeSystem() = 0;
	int threadWaitingLoop(int thrIdx, t_threadCode code, CThreadEnumerator *threadEnum, int nThread) const;
	virtual CRowSolution *FindSolution(size_t n, PERMUT_ELEMENT_TYPE lastRightPartIndex = -1)
															{ return NULL; }
	virtual bool TestFeatures(CEnumInfo *pEnumInfo)		    { return true; }
    virtual void prepareToTestExtraFeatures()               {}
	virtual void copyInfoFromMaster(const CEnumerator *pMaster) {}
    virtual CColOrbit **getUnforceColOrbPntr() const        { return NULL; }
#if USE_STRONG_CANONICITY_A
	void checkUnusedSolutions(CRowSolution *pRowSolution);
#else
	#define checkUnusedSolutions(x)
#endif
	void InitRowSolutions(const CEnumerator *pMaster);
	inline void setEnumInfo(CEnumInfo *pEnumInfo)			{ m_pEnumInfo = pEnumInfo; }
	inline CEnumInfo *enumInfo() const						{ return m_pEnumInfo; }
	inline bool isMatrOwner() const							{ return m_bMatrOwner; }

	const CMatrix *m_pMatrix;
	CRowSolution **m_pRow;
	CRowEquation *m_pRowEquation;
    CCanonicityChecker *m_pCanonChecker;
	bool m_bUseCanogGroup;
	CEnumInfo *m_pEnumInfo;
	bool m_bMatrOwner;
};
