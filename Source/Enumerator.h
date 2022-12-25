#pragma once
#include "RowSolution.h"
#include "InSysSolver.h"
#include "MatrixCanonChecker.h"
#include "DesignDB.h"

#define CONSTRUCTED_IN " constructed in"

#ifdef USE_CUDA
// Define this to turn on error checking
//#define CUDA_ERROR_CHECK
#define TRACE_ALL_CUDA_CALLS	0

#if !CONSTR_ON_GPU
#define releaseGPU_CanonChecker()			delete CanonCheckerGPU()
#define LAUNCH_CANONICITY_TESTING(x, y)		x->LaunchCanonicityTesting(y);
#else
#define setGPU_CanonChecker(x)
#define releaseGPU_CanonChecker()
#define LAUNCH_CANONICITY_TESTING(x, y)
#endif
#else
#define cudaDeviceReset()	// empty macro
#define cudaDeviceSetLimit(x,...)
#define setGPU_CanonChecker(x)
#define releaseGPU_CanonChecker()
#define LAUNCH_CANONICITY_TESTING(x, y)
#endif

#ifdef CUDA_ERROR_CHECK
#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

inline cudaError __cudaSafeCall(cudaError err, const char *file, const int line)
{
	if (TRACE_ALL_CUDA_CALLS || cudaSuccess != err)
	{
		fprintf(stderr, "\nSafeCall: %s at %s:%i\n",
			cudaSuccess != err? cudaGetErrorString(err) : "OK", file, line);
		if (cudaSuccess != err)
			exit(-1);
	}

	return err;
}

inline void __cudaCheckError(const char *file, const int line)
{
	cudaError err = cudaGetLastError();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "\nCheckErrorA %s at %s:%i\n",
			cudaSuccess != err ? cudaGetErrorString(err) : "OK", file, line);
		if (cudaSuccess != err)
			exit(-1);
	}
/*
	// More careful checking. However, this will affect performance.
	// Comment away if needed.
	err = cudaDeviceSynchronize();
	if (cudaSuccess != err)
	{
		fprintf(stderr, "\nCheckErrorB %s at %s:%i\n",
			cudaSuccess != err ? cudaGetErrorString(err) : "OK", file, line);
		if (cudaSuccess != err)
			exit(-1);
	} */
}
#else
#define CudaSafeCall( x )			x
#define CudaCheckError()
#endif

#define CudaSafeCallRet(x, y)		if (CudaSafeCall(x) != cudaSuccess) { return y; }

typedef enum {
	t_CPU,
	t_GPU
} t_ProcType;

Class2Def(CThreadEnumerator);

#if CONSTR_ON_GPU
#define MAKE_JOB_TITLE(x, y,...)
#else
#define MAKE_JOB_TITLE(x, y,...) x->makeJobTitle(y, __VA_ARGS__)
#endif

Class2Def(CEnumerator) : public Class2(CMatrixCanonChecker)
{
public:
	CK CEnumerator(const MatrixPntr pMatrix, uint enumFlags, int treadIdx = -1, uint nCanonChecker = 0);
	CC virtual ~CEnumerator();
	CK inline RowSolutionPntr rowStuff(T nRow = 0, T iPart = 0) const	{ return m_pRow[nRow] + iPart; }
	CK bool Enumerate(designParam *pParam, bool writeFile = false, EnumInfoPntr pEnumInfo = NULL, EnumeratorPntr pMaster = NULL, t_threadCode *pTreadCode = NULL);
	virtual bool makeJobTitle(const designParam *pParam, char *buffer, int len, const char *comment = "") const
															{ return false; }
	CK virtual T getX0_3() const							{ return 0; }
	CK inline RowSolutionPntr *rowStuffPntr() const			{ return m_pRow;  }
	CK virtual void setFirstUnforcedRow(T row = 0)			{}
	CK virtual T *forcibleLambdaPntr(T nRow = 0) const		{ return NULL; }
	CK virtual bool noReplicatedBlocks() const				{ return false; }
	CK virtual void CloneMasterInfo(const EnumeratorPntr p, size_t nRow) {}
	CK inline auto *designParams() const					{ return m_pParam; }
	CK inline T numRow() const								{ return nRow; }
	CK void setDesignDB(CDesignDB* pntr)					{ m_pDesignDB = pntr; }
	auto* designDB() const									{ return m_pDesignDB; }
#if CANON_ON_GPU
	CK inline auto CanonCheckerGPU() const					{ return m_pGPU_CanonChecker; }
	size_t copyColOrbitInfo(S nRow) const;
#endif
protected:
	CK virtual bool prepareToFindRowSolution()				{ return true; }
	CK inline InSysPntr getInSys() const					{ return this->IS_enumerator()? (InSysPntr)(this->matrix()) : NULL; }
	CK virtual void setX0_3(T value)						{}
	CK inline CSimpleArray<S>* rowEquation(T idx = 0) const	{ return m_pRowEquation + idx; }
	CK virtual bool checkSolutions(RowSolutionPntr p, S nPart, PERMUT_ELEMENT_TYPE idx, bool doSorting = true) { return false;  /* not implemented */ }
	CK inline void setRowEquation(CSimpleArray<S> *pntr)    { m_pRowEquation = pntr; }
	CK inline auto rowNumb() const							{ return this->matrix()->rowNumb(); }
	CK inline void setCurrentNumPart(T val)					{ m_nCurrentNumPart = val; }
	CK inline auto currentNumPart() const					{ return m_nCurrentNumPart; }
#if !CONSTR_ON_GPU
	virtual bool makeFileName(char *buffer, size_t len, const char *ext = NULL) const	{ return false; }
	bool getMasterFileName(char *buffer, size_t lenBuffer, size_t *pLenName) const;
#else
	#define makeFileName(buffer, ...)						false
#endif
	size_t getDirectory(char *buffer, size_t len, bool rowNeeded = true) const;
	CK inline void setUseCanonGroup(bool val)				{ m_bUseCanogGroup = val; }
	CK inline bool useCanonGroup() const					{ return m_bUseCanogGroup; }
	virtual void reset(T nRow, bool resetSolutions = true);
	CK void MakeRow(RowSolutionPntr pRowSolution, bool flag, S iFirstPartIdx = 0);
	CK virtual void CreateForcedRows()						{ this->setCurrentRowNumb(0); }
	CK virtual T firtstNonfixedRowNumber() const			{ return 2; }
	CK virtual bool fileExists(const char *path, bool file = true) const;
	CK virtual bool createNewFile(const char *fName) const	{ return true; }
	CK virtual bool SeekLogFile() const						{ return false; }

	CK virtual MatrixDataPntr CreateSpareMatrix(const EnumeratorPntr pMaster) { return NULL; }
	CK virtual void CreateAuxiliaryStructures(EnumeratorPntr pMaster)			{}
	CK virtual void InitGroupOderStorage(const CGroupOnParts<T>* pGroupOnParts)		{}
	CK virtual void ConstructedDesignProcessing() const		{}
	CK virtual void beforeEnumInfoOutput() const			{}
	CK void outBlockTitle(const char* title = "Constructed Matrices", bool checkFirstMatr = true) const;
	CK virtual void setDesignParams(designParam* ptr)		{ m_pParam = ptr; }
#if USE_MUTEX
	CK void setMaster(CEnumerator* pntr)					{ m_pMaster = pntr; }
	CK auto master() const									{ return m_pMaster; }
#else
#define setMaster(x)
#endif
private:
	virtual bool compareResults(char *fileName, size_t lenFileName, bool *pBetterResults = NULL);
	virtual void getEnumerationObjectKey(char *pKey, int len) const { strcpy_s(pKey, len, "EMPTY_KEY"); }
	virtual char *getEnumerationObjectKeyA(char* pKey, int len, const char* pKeyIn = NULL) { return NULL; }
	virtual int compareEnumerationDB_record(const char* record) { return 0; }
	virtual void outputTitle(FILE* file) const;
	virtual const char* getObjNameFormat() const			{ return "%9s:        "; }
	void UpdateEnumerationDB(char **pInfo, int len);
	bool cmpProcedure(FILE* file[2], bool* pBetterResults = NULL);
	CK virtual bool TestFeatures(EnumInfoPntr pEnumInfo, const MatrixDataPntr pMatrix, int *pMatrFlags = NULL, const EnumeratorPntr pEnum = NULL) const { return true; }
	CK virtual RowSolutionPntr setFirstRowSolutions()		{ return NULL; }
	CK RowSolutionPntr FindRowSolution(S *pPartNumb);
	CK virtual T MakeSystem(T numPart) = 0;
#if USE_THREADS
	int threadWaitingLoop(int thrIdx, t_threadCode code, ThreadEnumeratorPntr threadEnum, size_t nThread) const;
#endif
	CK virtual RowSolutionPntr FindSolution(T nVar, T nPart, PERMUT_ELEMENT_TYPE lastRightPartIndex = PERMUT_ELEMENT_MAX)
															{ return NULL; }
	CK virtual void prepareToTestExtraFeatures()			{}
	CK virtual void copyInfoFromMaster(const CEnumerator *pMaster) {}
	CK virtual int unforcedElement(const CColOrbit<S>* p, int nRow) const	{ return -1; }
	CK virtual ColOrbPntr* unforcedOrbits(T n, T idxPart) const	{ return NULL; }
	CK virtual void resetFirstUnforcedRow()					{}
	virtual S forcibleLambda(T nRow, T nPart) const			{ return ELEMENT_MAX; }
	virtual const char* getTopLevelDirName() const          { return NULL; }
	CK virtual void setFirstPartSolutionIndex(PERMUT_ELEMENT_TYPE idx) {}
	CK virtual PERMUT_ELEMENT_TYPE firstPartSolutionIndex(T nRow) const { return 0; }
	CK inline uchar *getSolutionsWereConstructed(T nParts, T rowNumb) const {
		return nParts > 1  && rowNumb < matrix()->rowNumb()? m_bSolutionsWereConstructed + rowNumb * nParts : NULL; }
	CK virtual void setForcibleLambda(T nRow, T val, T nPart) {}

#if PRINT_SOLUTIONS
	void printSolutions(const RowSolutionPntr pRowSolution, FILE* file, T nRow, bool markNextUsed, T nPartStart, T nPartEnd) const;
#endif

#if USE_STRONG_CANONICITY_A
	void checkUnusedSolutions(CRowSolution<T> *pRowSolution);
#else
	#define checkUnusedSolutions(x)
#endif
	CK void InitRowSolutions(const CEnumerator *pMaster);
#if CANON_ON_GPU
	CK inline void setGPU_CanonChecker(GPU_CanonChecker *pntr) { m_pGPU_CanonChecker = pntr; }
	bool TestCanonicityOnGPU() {
		int matrFlags = 0;
		if (!TestFeatures(enumInfo(), matrix(), &matrFlags))
			return true;

		assert(CanonCheckerGPU());

		return CanonCheckerGPU()->assignChecker(enumInfo(), this, matrFlags);
	}
#else
	#define TestCanonicityOnGPU()		false
#endif

	T nRow;
	RowSolutionPntr *m_pRow;
	CSimpleArray<S> *m_pRowEquation;
	bool m_bUseCanogGroup;
	designParam *m_pParam;
	ColOrbPntr* m_pFirstColOrb;
	T m_nCurrentNumPart;
	PERMUT_ELEMENT_TYPE* m_lastRightPartIndex;
	CDesignDB* m_pDesignDB = NULL;               // DB where the designs are stored

#if CANON_ON_GPU
	GPU_CanonChecker *m_pGPU_CanonChecker;
#endif
 protected:
	unsigned char* m_bSolutionsWereConstructed = NULL;
#if USE_MUTEX
	CEnumerator* m_pMaster = NULL;
	static std::mutex m_mutexDB;				 // mutex for accessing the database, it is used when BIBD  or the "master" for CBIBD is added
#endif
};

FClass2(CEnumerator)::CEnumerator(const MatrixPntr pMatrix, uint enumFlags, int treadIdx, uint nCanonChecker) :
	Class2(CMatrixCanonChecker)(pMatrix, enumFlags) {
	const auto numParts = this->numParts();
	m_pRow = new RowSolutionPntr[numParts * pMatrix->rowNumb()];
	m_lastRightPartIndex = new PERMUT_ELEMENT_TYPE[numParts];
	m_pFirstColOrb = new ColOrbPntr[numParts];
	setRowEquation(NULL);
	setGroupOnParts(NULL);
	setGPU_CanonChecker(nCanonChecker ? new Class2(CGPU_CanonChecker)(nCanonChecker, pMatrix, treadIdx) : NULL);
}

FClass2(CEnumerator)::~CEnumerator() {
#if USE_MUTEX
	const auto flag = master() == NULL;
	if (master()) {
		if (designParams()->thread_master_DB) {
			auto* pMasterDB = master()->designDB();
			m_mutexDB.lock();
			if (pMasterDB) {
				// Merging two Design DBs as ordered sets
				if (designDB()->recNumb()) {
					auto* pNewDB = new CDesignDB(pMasterDB->recordLength());
					pNewDB->combineDesignDBs(pMasterDB, designDB());
					master()->setDesignDB(pNewDB);
					delete pMasterDB;
				}
			}
			else {
				master()->setDesignDB(designDB());
				setDesignDB(NULL);
			}

			m_mutexDB.unlock();
			delete designDB();
		}
	}
	else {
		delete designDB();
	}
#endif

	delete[] rowStuff(this->rowMaster());
	delete[] rowStuffPntr();
	delete rowEquation();

	auto* pGroupOnParts = getGroupOnParts();
	if (pGroupOnParts && pGroupOnParts->owner() == this) {
		delete pGroupOnParts;
		setGroupOnParts(NULL);
	}

	delete[] m_lastRightPartIndex;
	delete[] m_pFirstColOrb;
	releaseGPU_CanonChecker();
}

template<>void CEnumerator<uchar,uchar>::rowSetFragm(uchar* pRow, uchar val, size_t len) const {
	memset(pRow, val, len);
}
