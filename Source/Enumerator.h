#pragma once
#include "RowSolution.h"
#include "InSysSolver.h"
#include "ColOrbitManager.h"
#include "matrix.h"

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

Class2Def(CMatrixCol) : public CColOrbitManager<S>
{
public:
	CC CMatrixCol(const MatrixDataPntr pMatrix, uint enumFlags = t_enumDefault) :
		CColOrbitManager<S>(pMatrix->maxElement() + 1, pMatrix->rowNumb(), pMatrix->colNumb(), pMatrix->numParts()) {
		initiateMatrixCol(pMatrix, enumFlags);
	}
	CC CMatrixCol(const MatrixDataPntr pMatrix, S rowNumb, S colNumb, S maxElem, uint enumFlags) :
		CColOrbitManager<S>(maxElem + 1, rowNumb, colNumb, pMatrix->numParts()) {
		initiateMatrixCol(pMatrix, enumFlags & (t_allFlags ^ t_matrixOwner));
	}
	CC ~CMatrixCol() {
		if (isMatrOwner())
			delete matrix();
	}
	CC inline void initiateMatrixCol(const MatrixDataPntr pMatrix, uint enumFlags = t_IS_enumerator) {
		m_pMatrix = pMatrix;
		setMatrOwner(enumFlags & t_matrixOwner);
		setIS_Enumerator(enumFlags & t_IS_enumerator);
		setOutFile(NULL);
	}
	CC inline const MatrixDataPntr matrix() const		{ return m_pMatrix; }
	CK inline bool IS_enumerator() const				{ return m_bIS_Emunerator; }
	inline void closeFile()								{ if (outFile()) { fclose(outFile()); setOutFile(NULL); } }
#if !CONSTR_ON_GPU
	CC inline void setOutFile(FILE *file)				{ m_pFile = file; }
#endif
	inline FILE *outFile() const						{ return m_pFile; }
	inline FILE **outFilePntr()							{ return &m_pFile; }
protected:
	CC inline void setIS_Enumerator(bool val)			{ m_bIS_Emunerator = val; }
private:
	CC inline void setMatrOwner(bool val)				{ m_bMatrOwner = val; }
	CC inline bool isMatrOwner() const					{ return m_bMatrOwner; }

	const MatrixDataPntr m_pMatrix;
	bool m_bMatrOwner;
	bool m_bIS_Emunerator;
	FILE *m_pFile;
};

#include "EnumInfo.h"

Class2Def(CMatrixCanonChecker) : public Class2(CMatrixCol), public Class2(CCanonicityChecker)
{
public:
	CC CMatrixCanonChecker(const MatrixDataPntr pMatrix, uint enumFlags) :
		Class2(CMatrixCol)(pMatrix, enumFlags),
		Class2(CCanonicityChecker)(pMatrix->rowNumb(), pMatrix->colNumb(), pMatrix->maxElement() + 1, enumFlags, pMatrix->numParts())
															{ setEnumInfo(NULL); }

	CC CMatrixCanonChecker(MatrixDataPntr pMatrix, S rowNumb, S colNumb, T maxElem, bool IS_enum) :
		Class2(CMatrixCol)(pMatrix, rowNumb, colNumb, maxElem, IS_enum),
		Class2(CCanonicityChecker)(rowNumb, colNumb, maxElem)	{ setEnumInfo(NULL); }
	CC ~CMatrixCanonChecker()								{ delete enumInfo(); }

	CC inline void setEnumInfo(EnumInfoPntr pEnumInfo)		{ m_pEnumInfo = pEnumInfo; }
	CC inline EnumInfoPntr enumInfo() const					{ return m_pEnumInfo; }
private:
	EnumInfoPntr m_pEnumInfo;
};

#if CANON_ON_GPU
Class2Def(CMatrixCanonCheckerGPU) : public Class2(CMatrixCanonChecker)
{
public:
	CC CMatrixCanonCheckerGPU(const Class2(CMatrixData) *pMatrix, uint enumFlags = 0) :
		Class2(CMatrixCanonChecker)(pMatrix, enumFlags)		{}
	CC ~CMatrixCanonCheckerGPU()							{ closeColOrbits(); }
};

#define TRACE_CUDA_FLAG		0
#if TRACE_CUDA_FLAG
#define FOPENN(x, y)		char buf[32]; sprintf_s(buf, "aaa_%d.txt", myID); FOPEN(file, buf, y); 
#define TRACE_CUDA(...)		{FOPENN(file, "a"); if (file) { \
								fprintf(file, "%f :", CGPU_CheckerInfo<MATRIX_ELEMENT_TYPE>::GetCurrentTime()); \
								fprintf(file, __VA_ARGS__); fclose(file); }}
#else
#define TRACE_CUDA(...)
#endif


typedef struct {
	COrderInfo* pOrderInfo;
	int nOrders;
	int nOrdersMax;
} GroupOrderInfo;

TFunc2(MakeCopyGroupInfo, void) (MatrixCanonCheckerGPU** ppCheckers, GroupOrderInfo* orderInfo,
	int CPU_threadIdx, cudaStream_t stream, COrderInfo* pOrderInfo);

Class2Def(CEnumerator);

TFunc2(AssignChecker, bool)(MatrixCanonCheckerGPU** ppCheckers, Class2(CMatrixData)** pAllMatrData, uint checkerIdx, const Enumerator* pCanonChecker, cudaStream_t stream
#if TRACE_CUDA_FLAG
	, int myID
#endif
);

TFunc2(ResetEnumInfoGlobal, void)(MatrixCanonCheckerGPU** ppCheckers, cudaStream_t stream);

Class2Def(CGPU_CheckerInfo)
{
public:
	static void InitCanonInfo(size_t nCPUthreads) {
		CudaSafeCall(cudaMalloc(&m_pOrderInfo, (m_nCPUthreads = nCPUthreads) * sizeof(m_pOrderInfo[0])));
		m_ppOrderInfo = new COrderInfo *[nCPUthreads];
		ResetCanonInfo();
		m_timer.startClock();
	}

	static void ResetCanonInfo() {
		CudaSafeCall(cudaMemset(m_pOrderInfo, 0, m_nCPUthreads * sizeof(m_pOrderInfo[0])));
		memset(m_ppOrderInfo, 0, m_nCPUthreads * sizeof(m_ppOrderInfo[0]));
	}

	static void CloseCanonInfo() {
		for (auto i = m_nCPUthreads; i--;)
			CudaSafeCall(cudaFree(m_ppOrderInfo[i]));

		delete[] m_ppOrderInfo;
		CudaSafeCall(cudaFree(m_pOrderInfo)); 
	}

	static float GetCurrentTime() { m_timer.setRunTime();  return m_timer.runTime(); }
	static COrderInfo *CopyCanonInfo(MatrixCanonCheckerGPU **ppCheckers, int threadID, int *nOrders, cudaStream_t stream
#if TRACE_CUDA_FLAG
		, int myID
#endif
	) {
		GroupOrderInfo groupInfo;
		COrderInfo *pOrderInfo = NULL;
		while (true) {
			MakeCopyGroupInfo(ppCheckers, getGroupOrderInfo(), threadID, stream, pOrderInfo);
			TRACE_CUDA(" CopyCanonInfo ==> threadID = %d  m_pOrderInfo = %p\n", threadID, m_pOrderInfo);
			CudaSafeCall(cudaMemcpyAsync(&groupInfo, m_pOrderInfo + threadID, sizeof(groupInfo), cudaMemcpyDeviceToHost, stream));

			*nOrders = groupInfo.nOrders;
			if (!*nOrders)
				return NULL;

			if (*nOrders > 0)
				break;

			CudaSafeCall(cudaFree(groupInfo.pOrderInfo));
			CudaSafeCall(cudaMalloc(&pOrderInfo, -*nOrders * sizeof(COrderInfo)));
			m_ppOrderInfo[threadID] = pOrderInfo;
		}

		COrderInfo *pOrderInfoCPU = new COrderInfo[*nOrders];
		TRACE_CUDA(" CopyCanonInfo ==> threadID = %d  (*nOrders) = %d\n", threadID, *nOrders);
		CudaSafeCall(cudaMemcpyAsync(pOrderInfoCPU, groupInfo.pOrderInfo, *nOrders * sizeof(pOrderInfo[0]), cudaMemcpyDeviceToHost, stream));
		return pOrderInfoCPU;
	}
	static GroupOrderInfo *getGroupOrderInfo() { return m_pOrderInfo; }
private:
	static size_t m_nCPUthreads;
	static GroupOrderInfo *m_pOrderInfo;
	static COrderInfo **m_ppOrderInfo;
	static CTimerInfo m_timer;
};

#define InitCanonInfo(x)	CGPU_CheckerInfo<TDATA_TYPES>::InitCanonInfo(x)
#define CloseCanonInfo()	CGPU_CheckerInfo<TDATA_TYPES>::CloseCanonInfo(); cudaDeviceReset();

Class2Def(CGPU_CanonChecker) : public GPU_CheckerInfo
{
public:
	CGPU_CanonChecker(uint total, const Class2(CMatrixData) *pMatrix, int treadIdx) : m_CPUtreadIdx(treadIdx), m_nCheckersTotal(total) {
		m_nCheckersAssigned = 0;
		m_pMatrix = pMatrix;
		m_pCanonFlag = NULL;
		m_ppAllMatrData = NULL;
		m_pMatrFlagsGPU = m_pMatrFlags = NULL;
		m_pColOrbitData[t_CPU] = m_pColOrbitData[t_GPU] = NULL;
		m_ppCanonChecker = NULL;
#if TRACE_CUDA_FLAG
		extern int cntr;
		myID = cntr++;
		TRACE_CUDA(" CGPU_CanonChecker ==> %p  m_pMatrix = %p  colNumb = %d\n", this, m_pMatrix, m_pMatrix ? m_pMatrix->colNumb() : -999)
#endif
	}

	~CGPU_CanonChecker() {
		TRACE_CUDA("~CGPU_CanonChecker <== %p (%p) numCheckers() = %d  m_ppCanonChecker = %p  m_ppAllMatrData = %p\n", this, m_pMatrix, numCheckers(), m_ppCanonChecker, m_ppAllMatrData)
		if (m_pMatrix)
			return;

		void ReleaseCheckers(MatrixCanonCheckerGPU **ppCheckers, uint numCheckers, cudaStream_t stream);
		ReleaseCheckers(m_ppCanonChecker, numCheckers(), CUDA_stream());

		for (uint i = 0; i < numCheckers(); ++i) {
			if (!m_ppAllMatrData[i])
				break;

			CudaSafeCall(cudaFree(m_ppAllMatrData[i]));
		}

		delete[] m_ppAllMatrData;

		CudaSafeCall(cudaFree((void *)m_pCanonFlag));

		delete[] m_pMatrFlags;
		cudaFree(m_pMatrFlagsGPU);
		delete[] ColOrbitData(t_CPU);
		CudaSafeCall(cudaFree((void *)ColOrbitData(t_GPU)));
		cudaStreamDestroy(CUDA_stream());
	}

	bool assignChecker(Class2(CEnumInfo) *pEnumInfo, const Enumerator *pEnum, int matrFlags) {
		if (!assignChecker(pEnum, matrFlags))
			return false;

		if (startTesting())
			LaunchCanonicityTesting(pEnumInfo, pEnum);

		return true;
	}

	size_t *ColOrbitData(t_ProcType idx = t_CPU) const	{ return m_pColOrbitData[idx]; }
	inline cudaStream_t CUDA_stream() const				{ return m_CUDA_stream; }
	virtual bool noReplicatedBlocks() const				{ return false; }
	bool assignChecker(const Enumerator *pCanonChecker, int matrFlag) {
		TRACE_CUDA("  assignChecker_0 === this = %p  m_pMatrix = %p\n", this, m_pMatrix)
		if (m_pMatrix) {	// CUDA memory was not allocated yet
			cudaStreamCreate(&m_CUDA_stream);
			size_t len = numCheckers() * (sizeof(m_ppCanonChecker[0]) + sizeof(m_pCanonFlag[0]));
			void *pntr;
			const cudaError err = cudaMalloc(&pntr, len + sizeof(m_pGroupInfo[0]) * numCheckers());

			TRACE_CUDA("  assignChecker === this = %p  err = %d (%s)  len = %zu\n", this, err, cudaGetErrorString(err), len)
			CudaSafeCallRet(err, false);

			cudaMemset(m_pCanonFlag = (uchar *)pntr, 0, len);

			m_ppAllMatrData = new void *[numCheckers()];
			m_pMatrFlags = new int[numCheckers()];
			memset(m_ppAllMatrData, 0, numCheckers() * sizeof(m_ppAllMatrData[0]));

			CudaSafeCall(cudaMalloc(&m_pMatrFlagsGPU, numCheckers() * sizeof(m_pMatrFlagsGPU[0])));

			// For each ColOrbit we will need to keep no more than 2 numbers:
			m_pColOrbitData[t_CPU] = new size_t[2 * m_pMatrix->colNumb() * m_pMatrix->rowNumb()];

			m_ppCanonChecker = (MatrixCanonCheckerGPU **)(m_pCanonFlag + sizeof(m_pCanonFlag[0]) * numCheckers());
			m_pGroupInfo = (uint *)((char *)m_ppCanonChecker + numCheckers() * sizeof(m_ppCanonChecker[0]));
			CudaSafeCallRet(cudaMalloc(m_pColOrbitData + t_GPU, len * sizeof(m_pColOrbitData[0][0])), false);
			m_pMatrix = NULL;
		}

		// Allocate GPU memory for matrix, if needed
		TRACE_CUDA("  assignChecker 2: this = %p  m_ppAllMatrData[m_nCheckersAssigned = %d] %p\n", this, m_nCheckersAssigned, m_ppAllMatrData[m_nCheckersAssigned])
			if (!m_ppAllMatrData[m_nCheckersAssigned])
				CudaSafeCall(cudaMalloc(m_ppAllMatrData + m_nCheckersAssigned, pCanonChecker->matrix()->lenData() + sizeof(Class2(CMatrixData))));

		m_pMatrFlags[m_nCheckersAssigned] = matrFlag;

#if TRACE_CUDA_FLAG
		return AssignChecker(m_ppCanonChecker, (Class2(CMatrixData)**)m_ppAllMatrData, m_nCheckersAssigned++, pCanonChecker, CUDA_stream(), myID);
#else
		return AssignChecker(m_ppCanonChecker, (Class2(CMatrixData)**)m_ppAllMatrData, m_nCheckersAssigned++, pCanonChecker, CUDA_stream());
#endif
	}

	void LaunchCanonicityTesting(EnumInfo *pEnumInfo, const Enumerator *pEnum) const
	{
		TRACE_CUDA("  TestCanonicity: this = %p  m_nCheckersAssigned %d\n", this, m_nCheckersAssigned)
			if (!m_nCheckersAssigned)
				return;

		CudaSafeCall(cudaMemcpyAsync(m_pMatrFlagsGPU, m_pMatrFlags, m_nCheckersAssigned * sizeof(m_pMatrFlags[0]), 
					cudaMemcpyHostToDevice, CUDA_stream()));

		void TestCanonicity(uint nMatr, MatrixCanonCheckerGPU **ppCheckers, uchar *pCanonFlag, uint *pGroupInfo, int * pMatrFlag, cudaStream_t stream);
		TestCanonicity(m_nCheckersAssigned, m_ppCanonChecker, m_pCanonFlag, m_pGroupInfo, m_pMatrFlagsGPU, CUDA_stream());

		int nOrders;
		const COrderInfo *pOrderInfoBase = CopyCanonInfo(m_ppCanonChecker, m_CPUtreadIdx, &nOrders, CUDA_stream()
#if TRACE_CUDA_FLAG
			, myID
#endif
		);

		if (pOrderInfoBase) {
			pEnumInfo->RecalcCountersByGroupOrders(pOrderInfoBase, nOrders);
			delete[] pOrderInfoBase;
			ResetEnumInfoGlobal(m_ppCanonChecker, CUDA_stream());
		}

		m_nCheckersAssigned = 0;
	}
protected:
	inline bool startTesting() const		{ return m_nCheckersAssigned == numCheckers(); }
private:
	CK inline auto numCheckers() const			{ return m_nCheckersTotal; }

	uint *m_pGroupInfo;
	MatrixCanonCheckerGPU **m_ppCanonChecker;
	const Class2(CMatrixData) *m_pMatrix;
	uchar *m_pCanonFlag;
	int *m_pMatrFlags;
	int *m_pMatrFlagsGPU;
	void **m_ppAllMatrData;
	size_t *m_pColOrbitData[2];
	cudaStream_t m_CUDA_stream;
	const int m_CPUtreadIdx;
	const uint m_nCheckersTotal;
	mutable uint m_nCheckersAssigned;
#if TRACE_CUDA_FLAG
	int myID;
#endif
};
#else
#define InitCanonInfo(x)
#define CloseCanonInfo()
#endif

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
	CK inline RowSolutionPntr rowStuff(S nRow = 0, S iPart = 0) const	{ return m_pRow[nRow] + iPart; }
	CK bool Enumerate(designParam *pParam, bool writeFile = false, EnumInfoPntr pEnumInfo = NULL, const EnumeratorPntr pMaster = NULL, t_threadCode *pTreadCode = NULL);
	virtual bool makeJobTitle(const designParam *pParam, char *buffer, int len, const char *comment = "") const
															{ return false; }
	CK virtual S getX0_3() const							{ return 0; }
	CK inline RowSolutionPntr *rowStuffPntr() const			{ return m_pRow;  }
	CK virtual S firstUnforcedRow() const					{ return 0; }
	CK virtual void setFirstUnforcedRow(S row)				{}
	CK virtual S *forcibleLambdaPntr(S nRow = 0) const		{ return NULL; }
	CK virtual bool noReplicatedBlocks() const				{ return false; }
	CK virtual void CloneMasterInfo(const EnumeratorPntr p, size_t nRow) {}
	CK inline auto designParams() const						{ return m_pParam; }

#if CANON_ON_GPU
	CK inline auto CanonCheckerGPU() const					{ return m_pGPU_CanonChecker; }
	size_t copyColOrbitInfo(S nRow) const;
#endif
protected:
	CK virtual bool prepareToFindRowSolution()				{ return true; }
	CK inline InSysPntr getInSys() const					{ return this->IS_enumerator()? (InSysPntr)(this->matrix()) : NULL; }
	CK virtual void setX0_3(S value)						{}
	CK inline CSimpleArray<S>* rowEquation(S idx = 0) const	{ return m_pRowEquation + idx; }
	CK virtual bool checkSolutions(RowSolutionPntr p, S nPart, PERMUT_ELEMENT_TYPE idx, bool doSorting = true) { return false;  /* not implemented */ }
	CK inline void setRowEquation(CSimpleArray<S> *pntr)    { m_pRowEquation = pntr; }
	CK inline S rowNumb() const								{ return this->matrix()->rowNumb(); }
	CK inline void setCurrentNumPart(S val)					{ m_nCurrentNumPart = val; }
	CK inline S currentNumPart() const						{ return m_nCurrentNumPart; }
#if !CONSTR_ON_GPU
	virtual bool makeFileName(char *buffer, size_t len, const char *ext = NULL) const	{ return false; }
	bool getMasterFileName(char *buffer, size_t lenBuffer, size_t *pLenName) const;
#else
	#define makeFileName(buffer, ...)						false
#endif
	size_t getDirectory(char *buffer, size_t len, bool rowNeeded = true) const;
	CK inline void setUseCanonGroup(bool val)				{ m_bUseCanogGroup = val; }
	CK inline bool useCanonGroup() const					{ return m_bUseCanogGroup; }
	virtual void reset(S nRow);
	CK ColOrbPntr MakeRow(const S *pRowSolution, bool nextColOrbNeeded, S partIdx = 0) const;
	CK void MakeRow(RowSolutionPntr pRowSolution, bool flag, S iFirstPartIdx = 0);
	CK virtual void CreateForcedRows()						{ this->setCurrentRowNumb(0); }
	CK virtual S firtstNonfixedRowNumber() const			{ return 2; }
	CK virtual bool fileExists(const char *path, bool file = true) const;
	CK virtual bool createNewFile(const char *fName) const	{ return true; }
	CK virtual bool SeekLogFile() const						{ return false; }
	void rowSetFragm(T *pRow, T val, size_t len) const {
		for (auto j = len; j--;)
			pRow[j] = val;
	}
private:
	virtual bool compareResults(char *fileName, size_t lenFileName, bool *pBetterResults = NULL) const;
	virtual void getEnumerationObjectKey(char *pInfo, int len) const { strcpy_s(pInfo, len, "EMPTY_KEY"); }
	virtual void outputTitle(FILE* file) const;
	virtual const char* getObjNameFormat() const			{ return "%9s:        "; }
	void UpdateEnumerationDB(char **pInfo, int len) const;
	bool cmpProcedure(FILE* file[2], bool* pBetterResults = NULL) const;
	CK virtual bool TestFeatures(EnumInfoPntr pEnumInfo, const MatrixDataPntr pMatrix, int *pMatrFlags = NULL, EnumeratorPntr pEnum = NULL) const { return true; }
	CK virtual RowSolutionPntr setFirstRowSolutions()		{ return NULL; }
	CK RowSolutionPntr FindRowSolution(S *pPartNumb);
	CK virtual S MakeSystem(S numPart) = 0;
#if USE_THREADS
	int threadWaitingLoop(int thrIdx, t_threadCode code, ThreadEnumeratorPntr threadEnum, size_t nThread) const;
#endif
	CK virtual RowSolutionPntr FindSolution(S nVar, S nPart, PERMUT_ELEMENT_TYPE lastRightPartIndex = PERMUT_ELEMENT_MAX)
															{ return NULL; }
	CK virtual void prepareToTestExtraFeatures()			{}
	CK virtual void copyInfoFromMaster(const CEnumerator *pMaster) {}
	CK virtual CColOrbit<S> **getUnforcedColOrbPntr(S idxPart) const	{ return NULL; }
	CK virtual int unforcedElement(const CColOrbit<S>* p, int nRow) const	{ return -1; }
	CK virtual ColOrbPntr* unforcedOrbits(S n, S idxPart) const	{ return NULL; }
	CK virtual void resetFirstUnforcedRow()					{}
	virtual S forcibleLambda(S nRow, S nPart) const			{ return ELEMENT_MAX; }
	virtual const char* getTopLevelDirName() const          { return NULL; }
	inline void setDesignParams(designParam* pntr)			{ m_pParam = pntr; }
#if PRINT_SOLUTIONS
	void printSolutions(const RowSolutionPntr pRowSolution, FILE* file, S nRow, bool markNextUsed, S nPart) const;
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

	RowSolutionPntr *m_pRow;
	CSimpleArray<S> *m_pRowEquation;
	bool m_bUseCanogGroup;
	designParam *m_pParam;
	ColOrbPntr* m_pFirstColOrb;
	S m_nCurrentNumPart;
	PERMUT_ELEMENT_TYPE* m_lastRightPartIndex;
#if CANON_ON_GPU
	GPU_CanonChecker *m_pGPU_CanonChecker;
#endif
};

FClass2(CEnumerator)::CEnumerator(const MatrixPntr pMatrix, uint enumFlags, int treadIdx, uint nCanonChecker) :
	Class2(CMatrixCanonChecker)(pMatrix, enumFlags) {
	const auto numParts = this->numParts();
	m_pRow = new RowSolutionPntr[numParts * pMatrix->rowNumb()];
	m_lastRightPartIndex = new PERMUT_ELEMENT_TYPE[numParts];
	m_pFirstColOrb = new ColOrbPntr[numParts];
	setRowEquation(NULL);
	setGPU_CanonChecker(nCanonChecker ? new Class2(CGPU_CanonChecker)(nCanonChecker, pMatrix, treadIdx) : NULL);
}

FClass2(CEnumerator)::~CEnumerator() {
	delete[] rowStuff(this->rowMaster());
	delete[] rowStuffPntr();
	if (false && numParts() > 1) // Will not use for now  (search for this comment)
		delete [] rowEquation();
	else
		delete rowEquation();

	delete[] m_lastRightPartIndex;
	delete[] m_pFirstColOrb;
	releaseGPU_CanonChecker();
}

template<>void CEnumerator<uchar,uchar>::rowSetFragm(uchar* pRow, uchar val, size_t len) const {
	memset(pRow, val, len);
}