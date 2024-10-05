#pragma once
#include "ColOrbitManager.h"
#include "matrix.h"
#include "EnumInfo.h"

Class2Def(CEnumInfo);

Class2Def(CMatrixCol) : public CColOrbitManager<S>
{
public:
	CC CMatrixCol(const InSysPntr pMatrix, uint enumFlags = t_enumDefault) :
		CColOrbitManager<S>(pMatrix->maxElement() + 1, pMatrix->rowNumb(), pMatrix->colNumb(), pMatrix->numParts()) {
		initiateMatrixCol(pMatrix, enumFlags);
	}
	CC CMatrixCol(const InSysPntr pMatrix, T rowNumb, T colNumb, T maxElem, uint enumFlags) :
		CColOrbitManager<S>(maxElem + 1, rowNumb, colNumb, pMatrix->numParts()) {
		initiateMatrixCol(pMatrix, enumFlags & (t_allFlags ^ t_matrixOwner));
	}
	CC ~CMatrixCol() {
		if (isMatrOwner())
			delete matrix();
	}
	CC inline void initiateMatrixCol(const InSysPntr pMatrix, uint enumFlags = t_IS_enumerator) {
		m_pMatrix = pMatrix;
		setMatrOwner(enumFlags & t_matrixOwner);
		setIS_Enumerator(enumFlags & t_IS_enumerator);
		setOutFile(NULL);
	}
	CC inline const auto matrix() const					{ return m_pMatrix; }
	CK inline bool IS_enumerator() const				{ return m_bIS_Emunerator; }
	inline void closeFile()								{ if (outFile()) { fclose(outFile()); setOutFile(NULL); } }
#if !CONSTR_ON_GPU
	CC inline void setOutFile(FILE * file)				{ m_pFile = file; }
#endif
	inline FILE * outFile() const						{ return m_pFile; }
	inline FILE * *outFilePntr()						{ return &m_pFile; }
protected:
	CC inline void setIS_Enumerator(bool val)			{ m_bIS_Emunerator = val; }
	void rowSetFragm(S* pRow, S val, size_t len) const {
		for (auto j = len; j--;)
			pRow[j] = val;
	}
	CK virtual CColOrbit<S>** getUnforcedColOrbPntr(T idxPart) const { return NULL; }
	CK virtual T firstUnforcedRow() const				{ return 0; }
private:
	CC inline void setMatrOwner(bool val)				{ m_bMatrOwner = val; }
	CC inline bool isMatrOwner() const					{ return m_bMatrOwner; }

	const InSysPntr m_pMatrix;
	bool m_bMatrOwner;
	bool m_bIS_Emunerator;
	FILE* m_pFile;
};

Class2Def(CMatrixCanonChecker) : public Class2(CMatrixCol), public Class2(CCanonicityChecker)
{
public:
	CC CMatrixCanonChecker(const InSysPntr pMatrix, uint enumFlags) :
		Class2(CMatrixCol)(pMatrix, enumFlags),
		Class2(CCanonicityChecker)(pMatrix->rowNumb(), pMatrix->colNumb(), pMatrix->maxElement() + 1, enumFlags, pMatrix->numParts())
	{}

	CC CMatrixCanonChecker(InSysPntr pMatrix, T rowNumb, T colNumb, S maxElem, bool IS_enum) :
		Class2(CMatrixCol)(pMatrix, rowNumb, colNumb, maxElem, IS_enum),
		Class2(CCanonicityChecker)(rowNumb, colNumb, maxElem)	{}
	CC ~CMatrixCanonChecker();

	CC inline void setEnumInfo(EnumInfoPntr pEnumInfo)			{ m_pEnumInfo = pEnumInfo; }
	CC inline EnumInfoPntr enumInfo() const						{ return m_pEnumInfo; }
	CK ColOrbPntr MakeRow(T nRow, const T* pRowSolution, uint clean_flags = t_MatrixFlags::t_getNextColOrb, T partIdx = 0) const;
	CK bool CheckBlockIntersections(T nRow, T b, const T* pRowSolution, T* pBlockIdx, T partIdx = 0);

	CK void GenerateBinaryColumnOrbits(T nRow, S *pRow, T *pColPermut = NULL) const;
	CK void GenerateColumnOrbits(T nRow, S* pRow, T* pColPermut = NULL) const;
	CK S* CanonizeMatrix(int k = 0, CanonicityCheckerPntr* ppClassGroup = NULL, T numClasses = 0);
	void sortRowsUpdateColumnOrbits(T v, T b, T nRowStart, bool initFlag = false);
	void adjustData(T rank, T colNumb) {
		CMatrixCol::setRank(rank);
		CCanonicityChecker::setRank(rank);
		setColNumber(colNumb);
	}
protected:
	CK inline auto* commonElemNumber() const							{ return m_pCommonElemNumber; }
	CK inline auto* blockIdx() const									{ return m_pBlockIdx; }
	CK inline auto* partIdx() const										{ return m_pPartIdx; }
	CK inline void setCommonElemNumber(uchar *pntr)						{ m_pCommonElemNumber = pntr; }
	CK inline void setBlockIdx(T *pntr)									{ m_pBlockIdx = pntr; }
	CK inline void setPartIdx(T *pntr)									{ m_pPartIdx = pntr; }
	CK inline void setClassSize(T nBlocks)                              { m_nClassSize = nBlocks; }
	CK inline T classSize() const										{ return m_nClassSize; }
	CK void ResetBlockIntersections(T nRow, T partIdx);
private:

	EnumInfoPntr m_pEnumInfo = NULL;
	uchar* m_pCommonElemNumber = NULL;			 // Number of common elements for 2 blocks
	T* m_pBlockIdx = NULL;                       // Block indices containing each element
	T* m_pPartIdx = NULL;                        // Part indices which to store the information for block intercestion updates 
	T m_nClassSize = 0;                          // Number of blocks in one parallel class of k-System
};

#if CANON_ON_GPU
Class2Def(CMatrixCanonCheckerGPU) : public Class2(CMatrixCanonChecker)
{
public:
	CC CMatrixCanonCheckerGPU(const Class2(CMatrixData) * pMatrix, uint enumFlags = 0) :
		Class2(CMatrixCanonChecker)(pMatrix, enumFlags) {}
	CC ~CMatrixCanonCheckerGPU() { closeColOrbits(); }
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
		m_ppOrderInfo = new COrderInfo * [nCPUthreads];
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
	static COrderInfo* CopyCanonInfo(MatrixCanonCheckerGPU * *ppCheckers, int threadID, int* nOrders, cudaStream_t stream
#if TRACE_CUDA_FLAG
		, int myID
#endif
	) {
		GroupOrderInfo groupInfo;
		COrderInfo* pOrderInfo = NULL;
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

		COrderInfo* pOrderInfoCPU = new COrderInfo[*nOrders];
		TRACE_CUDA(" CopyCanonInfo ==> threadID = %d  (*nOrders) = %d\n", threadID, *nOrders);
		CudaSafeCall(cudaMemcpyAsync(pOrderInfoCPU, groupInfo.pOrderInfo, *nOrders * sizeof(pOrderInfo[0]), cudaMemcpyDeviceToHost, stream));
		return pOrderInfoCPU;
	}
	static GroupOrderInfo* getGroupOrderInfo() { return m_pOrderInfo; }
private:
	static size_t m_nCPUthreads;
	static GroupOrderInfo* m_pOrderInfo;
	static COrderInfo** m_ppOrderInfo;
	static CTimerInfo m_timer;
};

#define InitCanonInfo(x)	CGPU_CheckerInfo<TDATA_TYPES>::InitCanonInfo(x)
#define CloseCanonInfo()	CGPU_CheckerInfo<TDATA_TYPES>::CloseCanonInfo(); cudaDeviceReset();

Class2Def(CGPU_CanonChecker) : public GPU_CheckerInfo
{
public:
	CGPU_CanonChecker(uint total, const Class2(CMatrixData) * pMatrix, int treadIdx) : m_CPUtreadIdx(treadIdx), m_nCheckersTotal(total) {
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

		void ReleaseCheckers(MatrixCanonCheckerGPU * *ppCheckers, uint numCheckers, cudaStream_t stream);
		ReleaseCheckers(m_ppCanonChecker, numCheckers(), CUDA_stream());

		for (uint i = 0; i < numCheckers(); ++i) {
			if (!m_ppAllMatrData[i])
				break;

			CudaSafeCall(cudaFree(m_ppAllMatrData[i]));
		}

		delete[] m_ppAllMatrData;

		CudaSafeCall(cudaFree((void*)m_pCanonFlag));

		delete[] m_pMatrFlags;
		cudaFree(m_pMatrFlagsGPU);
		delete[] ColOrbitData(t_CPU);
		CudaSafeCall(cudaFree((void*)ColOrbitData(t_GPU)));
		cudaStreamDestroy(CUDA_stream());
	}

	bool assignChecker(Class2(CEnumInfo) * pEnumInfo, const Enumerator * pEnum, int matrFlags) {
		if (!assignChecker(pEnum, matrFlags))
			return false;

		if (startTesting())
			LaunchCanonicityTesting(pEnumInfo, pEnum);

		return true;
	}

	size_t * ColOrbitData(t_ProcType idx = t_CPU) const { return m_pColOrbitData[idx]; }
	inline cudaStream_t CUDA_stream() const { return m_CUDA_stream; }
	virtual bool noReplicatedBlocks() const { return false; }
	bool assignChecker(const Enumerator * pCanonChecker, int matrFlag) {
		TRACE_CUDA("  assignChecker_0 === this = %p  m_pMatrix = %p\n", this, m_pMatrix)
		if (m_pMatrix) {	// CUDA memory was not allocated yet
			cudaStreamCreate(&m_CUDA_stream);
			size_t len = numCheckers() * (sizeof(m_ppCanonChecker[0]) + sizeof(m_pCanonFlag[0]));
			void* pntr;
			const cudaError err = cudaMalloc(&pntr, len + sizeof(m_pGroupInfo[0]) * numCheckers());

			TRACE_CUDA("  assignChecker === this = %p  err = %d (%s)  len = %zu\n", this, err, cudaGetErrorString(err), len)
			CudaSafeCallRet(err, false);

			cudaMemset(m_pCanonFlag = (uchar*)pntr, 0, len);

			m_ppAllMatrData = new void* [numCheckers()];
			m_pMatrFlags = new int[numCheckers()];
			memset(m_ppAllMatrData, 0, numCheckers() * sizeof(m_ppAllMatrData[0]));

			CudaSafeCall(cudaMalloc(&m_pMatrFlagsGPU, numCheckers() * sizeof(m_pMatrFlagsGPU[0])));

			// For each ColOrbit we will need to keep no more than 2 numbers:
			m_pColOrbitData[t_CPU] = new size_t[2 * m_pMatrix->colNumb() * m_pMatrix->rowNumb()];

			m_ppCanonChecker = (MatrixCanonCheckerGPU**)(m_pCanonFlag + sizeof(m_pCanonFlag[0]) * numCheckers());
			m_pGroupInfo = (uint*)((char*)m_ppCanonChecker + numCheckers() * sizeof(m_ppCanonChecker[0]));
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

	void LaunchCanonicityTesting(EnumInfo * pEnumInfo, const Enumerator * pEnum) const
	{
		TRACE_CUDA("  TestCanonicity: this = %p  m_nCheckersAssigned %d\n", this, m_nCheckersAssigned)
			if (!m_nCheckersAssigned)
				return;

		CudaSafeCall(cudaMemcpyAsync(m_pMatrFlagsGPU, m_pMatrFlags, m_nCheckersAssigned * sizeof(m_pMatrFlags[0]),
					cudaMemcpyHostToDevice, CUDA_stream()));

		void TestCanonicity(uint nMatr, MatrixCanonCheckerGPU * *ppCheckers, uchar * pCanonFlag, uint * pGroupInfo, int* pMatrFlag, cudaStream_t stream);
		TestCanonicity(m_nCheckersAssigned, m_ppCanonChecker, m_pCanonFlag, m_pGroupInfo, m_pMatrFlagsGPU, CUDA_stream());

		int nOrders;
		const COrderInfo* pOrderInfoBase = CopyCanonInfo(m_ppCanonChecker, m_CPUtreadIdx, &nOrders, CUDA_stream()
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
	inline bool startTesting() const { return m_nCheckersAssigned == numCheckers(); }
private:
	CK inline auto numCheckers() const { return m_nCheckersTotal; }

	uint * m_pGroupInfo;
	MatrixCanonCheckerGPU** m_ppCanonChecker;
	const Class2(CMatrixData)* m_pMatrix;
	uchar* m_pCanonFlag;
	int* m_pMatrFlags;
	int* m_pMatrFlagsGPU;
	void** m_ppAllMatrData;
	size_t* m_pColOrbitData[2];
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

