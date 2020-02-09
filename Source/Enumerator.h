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
#define releaseGPU_CanonChecker() delete GPU_CanonChecker()
#define LAUNCH_CANONICITY_TESTING(x, y) x->LaunchCanonicityTesting(y);
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

template<class T>
class CMatrixCol : public CColOrbitManager<T>
{
public:
	CC CMatrixCol(const CMatrixData<T> *pMatrix, uint enumFlags = t_enumDefault) :
		CColOrbitManager<T>(pMatrix->maxElement() + 1, pMatrix->rowNumb(), pMatrix->colNumb()) {
		initiateMatrixCol(pMatrix, enumFlags);
	}
	CC CMatrixCol(const CMatrixData<T> *pMatrix, T rowNumb, T colNumb, T maxElem, uint enumFlags) :
		CColOrbitManager<T>(maxElem + 1, rowNumb, colNumb) {
		initiateMatrixCol(pMatrix, enumFlags & (t_allFlags ^ t_matrixOwner));
	}
	CC ~CMatrixCol() {
		if (isMatrOwner())
			delete static_cast<const CMatrix<T> *>(matrix());
	}
	CC inline void initiateMatrixCol(const CMatrixData<T> *pMatrix, uint enumFlags = t_IS_enumerator) {
		m_pMatrix = pMatrix;
		setMatrOwner(enumFlags & t_matrixOwner);
		setIS_Enumerator(enumFlags & t_IS_enumerator);
		setOutFile(NULL);
	}
	CC inline const CMatrixData<T> *matrix() const	{ return m_pMatrix; }
	CK inline bool IS_enumerator() const			{ return m_bIS_Emunerator; }
	inline void closeFile()							{ if (outFile()) { fclose(outFile()); setOutFile(NULL); } }
#if !CONSTR_ON_GPU
	CC inline void setOutFile(FILE *file)			{ m_pFile = file; }
#endif
	inline FILE *outFile() const					{ return m_pFile; }
	inline FILE **outFilePntr()						{ return &m_pFile; }
protected:
	CC inline void setIS_Enumerator(bool val)		{ m_bIS_Emunerator = val; }
private:
	CC inline void setMatrOwner(bool val)			{ m_bMatrOwner = val; }
	CC inline bool isMatrOwner() const				{ return m_bMatrOwner; }

	const CMatrixData<T> *m_pMatrix;
	bool m_bMatrOwner;
	bool m_bIS_Emunerator;
	FILE *m_pFile;
};

#include "EnumInfo.h"

template<class T>
class CMatrixCanonChecker : public CMatrixCol<T>, public CCanonicityChecker<T>
{
public:
	CC CMatrixCanonChecker(const CMatrixData<T> *pMatrix, uint enumFlags) :
		CMatrixCol<T>(pMatrix, enumFlags),
		CCanonicityChecker<T>(pMatrix->rowNumb(), pMatrix->colNumb(), pMatrix->maxElement() + 1, enumFlags)
															{ setEnumInfo(NULL); }

	CC CMatrixCanonChecker(CMatrixData<T> *pMatrix, T rowNumb, T colNumb, T maxElem, bool IS_enum) :
		CMatrixCol<T>(pMatrix, rowNumb, colNumb, maxElem, IS_enum),
		CCanonicityChecker<T>(rowNumb, colNumb, maxElem)	{ setEnumInfo(NULL); }
	CC ~CMatrixCanonChecker()								{ delete enumInfo(); }

	CC inline void setEnumInfo(CEnumInfo<T> *pEnumInfo)		{ m_pEnumInfo = pEnumInfo; }
	CC inline CEnumInfo<T> *enumInfo() const				{ return m_pEnumInfo; }
private:
	CEnumInfo<T> *m_pEnumInfo;
};

#if CANON_ON_GPU
template<class T>
class CMatrixCanonCheckerGPU : public CMatrixCanonChecker<T>
{
public:
	CC CMatrixCanonCheckerGPU(const CMatrixData<T> *pMatrix, uint enumFlags = 0) :
		CMatrixCanonChecker<T>(pMatrix, enumFlags)			{}
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
	COrderInfo *pOrderInfo;
	int nOrders;
	int nOrdersMax;
} GroupOrderInfo;

template<class T>
class CGPU_CheckerInfo
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
	static COrderInfo *CopyCanonInfo(CMatrixCanonCheckerGPU<T> **ppCheckers, int threadID, int *nOrders, cudaStream_t stream
#if TRACE_CUDA_FLAG
		, int myID
#endif
	) {
		void MakeCopyGroupInfo(CMatrixCanonCheckerGPU<T> **ppCheckers, GroupOrderInfo *orderInfo, 
			int CPU_threadIdx, cudaStream_t stream, COrderInfo *pOrderInfo);
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

#define InitCanonInfo(x)	CGPU_CheckerInfo<MATRIX_ELEMENT_TYPE>::InitCanonInfo(x)
#define CloseCanonInfo()	CGPU_CheckerInfo<MATRIX_ELEMENT_TYPE>::CloseCanonInfo(); cudaDeviceReset();
template<class T> class CEnumerator;

template<class T>
class CGPU_CanonChecker : public CGPU_CheckerInfo<T>
{
public:
	CGPU_CanonChecker(uint total, const CMatrixData<T> *pMatrix, int treadIdx) : m_CPUtreadIdx(treadIdx), m_nCheckersTotal(total) {
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

		void ReleaseCheckers(CMatrixCanonCheckerGPU<T> **ppCheckers, uint numCheckers, cudaStream_t stream);
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

	bool assignChecker(CEnumInfo<T> *pEnumInfo, const CEnumerator<T> *pEnum, int matrFlags) {
		if (!assignChecker(pEnum, matrFlags))
			return false;

		if (startTesting())
			LaunchCanonicityTesting(pEnumInfo, pEnum);

		return true;
	}

	size_t *ColOrbitData(t_ProcType idx = t_CPU) const { return m_pColOrbitData[idx]; }
	inline cudaStream_t CUDA_stream() const { return m_CUDA_stream; }
	virtual bool noReplicatedBlocks() const { return false; }
	bool assignChecker(const CEnumerator<T> *pCanonChecker, int matrFlag) {
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

			m_ppCanonChecker = (CMatrixCanonCheckerGPU<T> **)(m_pCanonFlag + sizeof(m_pCanonFlag[0]) * numCheckers());
			m_pGroupInfo = (uint *)((char *)m_ppCanonChecker + numCheckers() * sizeof(m_ppCanonChecker[0]));
			CudaSafeCallRet(cudaMalloc(m_pColOrbitData + t_GPU, len * sizeof(m_pColOrbitData[0][0])), false);
			m_pMatrix = NULL;
		}

		// Allocate GPU memory for matrix, if needed
		TRACE_CUDA("  assignChecker 2: this = %p  m_ppAllMatrData[m_nCheckersAssigned = %d] %p\n", this, m_nCheckersAssigned, m_ppAllMatrData[m_nCheckersAssigned])
			if (!m_ppAllMatrData[m_nCheckersAssigned])
				CudaSafeCall(cudaMalloc(m_ppAllMatrData + m_nCheckersAssigned, pCanonChecker->matrix()->lenData() + sizeof(CMatrixData<T>)));

		m_pMatrFlags[m_nCheckersAssigned] = matrFlag;

		bool AssignChecker(CMatrixCanonCheckerGPU<T> **ppCheckers, CMatrixData<T> **pAllMatrData, uint checkerIdx, const CEnumerator<T> *pCanonChecker, cudaStream_t stream
#if TRACE_CUDA_FLAG
			, int myID
#endif
		);

#if TRACE_CUDA_FLAG
		return AssignChecker(m_ppCanonChecker, (CMatrixData<T> **)m_ppAllMatrData, m_nCheckersAssigned++, pCanonChecker, CUDA_stream(), myID);
#else
		return AssignChecker(m_ppCanonChecker, (CMatrixData<T> **)m_ppAllMatrData, m_nCheckersAssigned++, pCanonChecker, CUDA_stream());
#endif
	}

	void LaunchCanonicityTesting(CEnumInfo<T> *pEnumInfo, const CEnumerator<T> *pEnum) const
	{
		TRACE_CUDA("  TestCanonicity: this = %p  m_nCheckersAssigned %d\n", this, m_nCheckersAssigned)
			if (!m_nCheckersAssigned)
				return;

		CudaSafeCall(cudaMemcpyAsync(m_pMatrFlagsGPU, m_pMatrFlags, m_nCheckersAssigned * sizeof(m_pMatrFlags[0]), 
					cudaMemcpyHostToDevice, CUDA_stream()));

		void TestCanonicity(uint nMatr, CMatrixCanonCheckerGPU<T> **ppCheckers, uchar *pCanonFlag, uint *pGroupInfo, int * pMatrFlag, cudaStream_t stream);
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

			void ResetEnumInfoGlobal(CMatrixCanonCheckerGPU<T> **ppCheckers, cudaStream_t stream);
			ResetEnumInfoGlobal(m_ppCanonChecker, CUDA_stream());
		}

		m_nCheckersAssigned = 0;
	}
protected:
	inline bool startTesting() const		{ return m_nCheckersAssigned == numCheckers(); }
private:
	CK inline auto numCheckers() const			{ return m_nCheckersTotal; }

	uint *m_pGroupInfo;
	CMatrixCanonCheckerGPU<T> **m_ppCanonChecker;
	const CMatrixData<T> *m_pMatrix;
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

template<class T> class CThreadEnumerator;
#if CONSTR_ON_GPU
#define MAKE_JOB_TITLE(x, y,...)
#else
#define MAKE_JOB_TITLE(x, y,...) x->makeJobTitle(y, __VA_ARGS__)
#endif

template<class T>
class CEnumerator : public CMatrixCanonChecker<T>
{
public:
	CK CEnumerator(const CMatrix<T> *pMatrix, uint enumFlags, int treadIdx = -1, uint nCanonChecker = 0);
	CC virtual ~CEnumerator();
	CK inline CRowSolution<T> *rowStuff(size_t nRow = 0) const	{ return m_pRow[nRow]; }
	CK ulonglong Enumerate(designParam *pParam, bool writeFile = false, CEnumInfo<T> *pEnumInfo = NULL, const CEnumerator<T> *pMaster = NULL, t_threadCode *pTreadCode = NULL);
	virtual bool makeJobTitle(const designParam *pParam, char *buffer, int len, const char *comment = "") const
															{ return false; }
	CK virtual VECTOR_ELEMENT_TYPE getX0_3() const          { return 0; }
	CK inline CRowSolution<T> **rowStuffPntr() const						{ return m_pRow;  }
    CK virtual size_t firstUnforcedRow() const              { return 0; }
	CK virtual void setFirstUnforcedRow(size_t rowNum = 0)  {}
	CK virtual size_t *forcibleLambdaPntr() const           { return NULL; }
	CK virtual bool noReplicatedBlocks() const				{ return false; }
	CK virtual void CloneMasterInfo(const CEnumerator<T> *p, size_t nRow) {}
	CK inline designParam *designParams() const				{ return m_pParam; }

#if CANON_ON_GPU
	CK inline auto GPU_CanonChecker() const					{ return m_pGPU_CanonChecker; }
	size_t copyColOrbitInfo(T nRow) const;
#endif
protected:
	CK virtual bool prepareToFindRowSolution()				{ return true; }
	CK inline C_InSys<T> *getInSys() const					{ return this->IS_enumerator()? (C_InSys<T> *)this->matrix() : NULL; }
    CK virtual void setX0_3(VECTOR_ELEMENT_TYPE value)      {}
	CK inline CSimpleArray<T> *rowEquation() const			{ return m_pRowEquation; }
	virtual int unforcedElement(const CColOrbit<T> *p, int nRow) const    { return -1; }
    CK virtual bool sortSolutions(CRowSolution<T> *p, size_t idx) { return false;  /* not implemented */ }
	CK inline void setRowEquation(CSimpleArray<T> *pntr)    { m_pRowEquation = pntr; }
    CK inline T rowNumb() const								{ return this->matrix()->rowNumb(); }
#if !CONSTR_ON_GPU
	virtual bool makeFileName(char *buffer, size_t len, const char *ext = NULL) const	{ return false; }
	bool getMasterFileName(char *buffer, size_t lenBuffer, size_t *pLenName) const;
#else
	#define makeFileName(buffer, ...)						false
#endif
	size_t getDirectory(char *buffer, size_t len, bool rowNeeded = true) const;
	CK inline void setUseCanonGroup(bool val)				{ m_bUseCanogGroup = val; }
	CK inline bool useCanonGroup() const					{ return m_bUseCanogGroup; }
	virtual void reset(T nRow);
	CK CColOrbit<T> *MakeRow(const VECTOR_ELEMENT_TYPE *pRowSolution) const;
	CK virtual bool fileExists(const char *path, bool file = true) const;
	CK virtual bool createNewFile(const char *fName) const	{ return true; }
	CK virtual bool SeekLogFile() const						{ return false; }

private:
	virtual bool compareResults(char *fileName, size_t lenFileName, bool *pBetterResults = NULL) const;
	virtual void getEnumerationObjectKey(char *pInfo, int len) const { strcpy_s(pInfo, len, "EMPTY_KEY"); }
	virtual void outputTitle(FILE* file) const;
	virtual const char* getObjNameFormat() const			{ return "%9s:        "; }
	void UpdateEnumerationDB(char **pInfo, int len) const;
	bool cmpProcedure(FILE* file[2], bool* pBetterResults = NULL) const;
	CK virtual bool TestFeatures(CEnumInfo<T> *pEnumInfo, const CMatrixData<T> *pMatrix, int *pMatrFlags = NULL, CEnumerator<T> *pEnum = NULL) const { return true; }
	CK virtual CRowSolution<T> *setFirstRowSolutions()		{ return NULL; }
	CK CRowSolution<T> *FindRowSolution(PERMUT_ELEMENT_TYPE lastRightPartIndex = PERMUT_ELEMENT_MAX);
	CK virtual size_t MakeSystem() = 0;
#if USE_THREADS
	int threadWaitingLoop(int thrIdx, t_threadCode code, CThreadEnumerator<T> *threadEnum, size_t nThread) const;
#endif
	CK virtual CRowSolution<T> *FindSolution(size_t n, PERMUT_ELEMENT_TYPE lastRightPartIndex = PERMUT_ELEMENT_MAX)
															{ return NULL; }
    CK virtual void prepareToTestExtraFeatures()			{}
	CK virtual void copyInfoFromMaster(const CEnumerator *pMaster) {}
    CK virtual CColOrbit<T> **getUnforcedColOrbPntr() const	{ return NULL; }
	virtual size_t forcibleLambda(size_t i) const			{ return -1; }
	inline void setDesignParams(designParam *pntr)          { m_pParam = pntr; }
	virtual const char* getTopLevelDirName() const          { return NULL; }
#if USE_STRONG_CANONICITY_A
	void checkUnusedSolutions(CRowSolution<T> *pRowSolution);
#else
	#define checkUnusedSolutions(x)
#endif
	CK void InitRowSolutions(const CEnumerator *pMaster);
#if CANON_ON_GPU
	CK inline void setGPU_CanonChecker(CGPU_CanonChecker<T> *pntr) { m_pGPU_CanonChecker = pntr; }
	bool TestCanonicityOnGPU() {
		int matrFlags = 0;
		if (!TestFeatures(enumInfo(), matrix(), &matrFlags))
			return true;

		assert(GPU_CanonChecker());

		return GPU_CanonChecker()->assignChecker(enumInfo(), this, matrFlags);
	}
#else
	#define TestCanonicityOnGPU()		false
#endif

	CRowSolution<T> **m_pRow;
	CSimpleArray<T> *m_pRowEquation;
	bool m_bUseCanogGroup;
	designParam *m_pParam;
#if CANON_ON_GPU
	CGPU_CanonChecker<T> *m_pGPU_CanonChecker;
#endif
};

template<class T>
CEnumerator<T>::CEnumerator(const CMatrix<T> *pMatrix, uint enumFlags, int treadIdx, uint nCanonChecker) :
	CMatrixCanonChecker<T>(pMatrix, enumFlags)
{
	m_pRow = new CRowSolution<T> *[pMatrix->rowNumb()];
	setRowEquation(NULL);
	setGPU_CanonChecker(nCanonChecker ? new CGPU_CanonChecker<T>(nCanonChecker, pMatrix, treadIdx) : NULL);
}

template<class T>
CEnumerator<T>::~CEnumerator()
{
	delete[] rowStuff(this->rowMaster());
	delete[] rowStuffPntr();
	delete rowEquation();
	releaseGPU_CanonChecker();
}

