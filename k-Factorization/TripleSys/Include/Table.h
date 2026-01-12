#pragma once
#include "TripleSys.h"

template<typename T>
void outMatrix(const T* c, int nl, int nc, int np, int ns, FILE* f, bool makeString = false, 
	bool toScreen = false, const char *pStartLine = " \"", int cntr=-1, 
	const unsigned char* pDayPerm=NULL, bool empty_line = true) {
	char buffer[8192];
	const auto* endLine = makeString ? " \"\n" : "\n";
	const auto* elemFrmt = nc < 10 ? " %1d" : (nc < 100 ? " %2d" : " %3d");
	for (int j = 0; j < nl; j++) {
		char* pBuf = buffer;
		SPRINTFD(pBuf, buffer, pStartLine);
		for (int i = 0; i < nc; i++, c++) {
			if (np && !(i % np))
				SPRINTFD(pBuf, buffer, " ");
			/* need fix
			if (*c == -1 && !f)
				printfGreen(" %3d", *c);
			else */
				SPRINTFD(pBuf, buffer, elemFrmt, *c);
		}

		if (cntr < 0) {
			if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
				SPRINTFD(pBuf, buffer, endLine);
			else
				SPRINTFD(pBuf, buffer, " ");
		}
		else {
			if (cntr)
				SPRINTFD(pBuf, buffer, "\":  day =%2d\n", pDayPerm[j]);
			else
				SPRINTFD(pBuf, buffer, "\"\n");
		}

		_printf(f, toScreen, buffer);
	}

	if (empty_line)
		_printf(f, toScreen, "\n");
}

template<typename T>
class IOutGroupHandle {
public:
	virtual void makeGroupOutput(const CRepository<T>* pElemInfo, bool outToScreen = false, bool checkNestedGroups = true) = 0;
	virtual ~IOutGroupHandle() {}
	inline void setOutFileName(const char* pFileName, bool resetFile = true) {
		if ((m_pFileName = pFileName) && resetFile)
			std::remove(pFileName);
	}
protected:
	const auto outFileName() const				{ return m_pFileName; }
	const char* m_pFileName = NULL;
};

template<typename T>
class Table : public IOutGroupHandle<T> {
public:
	Table(char const *name, int nl, int nc, int ns, int np, bool makeString = false, bool outCntr = false) :
		m_pName(name), m_nc(nc), m_nl(nl), m_ns(ns), m_np(np), 
		m_makeString(makeString), m_bOutCntr(outCntr), IOutGroupHandle<T>() {}
	void printTable(const T *c, bool outCntr = false, bool outToScreen = true, int nl = 0, const char* pStartLine = " \"", const int* idx = NULL);
	void printTableInfo(const char* pInfo);
	inline void addCounterToTableName(bool val) { m_bOutCntr = val; }
	inline auto counter() const					{ return m_cntr; }
	virtual const char *name() 					{ return m_pName; }
	inline void setName(const char* pName)		{ m_pName = pName; }
	virtual void makeGroupOutput(const CRepository<T>* pElemInfo, bool outToScreen = false, bool checkNestedGroups = true) {
		this->printTable((const T*)pElemInfo->getObject(), false, outToScreen, pElemInfo->orderOfGroup(), "", pElemInfo->getIndices());
	}
	virtual void outLS(ctchar* c, int nl, FILE* f, bool outToScreen) {}
	virtual bool isLSCreated() const			{ return false; }
protected:
	inline auto groupSize() const				{ return m_np; }

	const char * m_pName;
	const ushort m_nc;
	const int m_nl;
	const int m_ns;
	const int m_np;
private:
	const bool m_makeString;
	bool m_bOutCntr;          // When true, the counter will be added to the Table Name
public:
	size_t m_cntr = 0;
};

class TableAut : public Table<tchar>
{
public:
	TableAut(char const* name, int nl, int nc, int ns, int np, bool makeString = false, bool outCntr = false) :
		Table<tchar>(name, nl, nc, ns, np, makeString, outCntr) {}
	virtual ~TableAut()							{ delete[] m_pBuffer; }
	inline void allocateBuffer(size_t len)		{ m_pBuffer = new char[m_nLenBuffer = len]; }
	inline void setGroupOrder(int groupOrder) 	{ m_groupOrder = groupOrder; }
	inline void setInfo(const char* pInfo)		{ m_pInfo = pInfo; }
	virtual const char* name() {
		const auto len = snprintf(m_pBuffer, m_nLenBuffer, "%s %d", m_pName, m_groupOrder);
		if (m_pInfo) {
			const auto new_len = len + strlen(m_pInfo) + 4;
			if (new_len > m_nLenBuffer) {
				const auto* pTmp = m_pBuffer;
				allocateBuffer(new_len);
				memcpy(m_pBuffer, pTmp, len);
				delete[] pTmp;
			}

			snprintf(m_pBuffer + len, new_len - len, ",  %s", m_pInfo);
			m_pInfo = NULL;
		}
		return m_pBuffer;
	}
private:
	int m_groupOrder = 0;
	char* m_pBuffer = NULL;
	const char* m_pInfo = NULL;
	size_t m_nLenBuffer = 0;
};

template<typename T>
class COutGroupHandle : public Table<T> {
public:
	COutGroupHandle(uint outGroupMask, char const* name, int degree) : m_outGroupMask(outGroupMask),
		Table<T>(name, 0, degree, -1, 0, false, true) {}
	virtual ~COutGroupHandle()	{}
protected: 
	const uint m_outGroupMask;  // In what forms and what the group's information should be printed
}; 

template<typename T>
class Generators : public CGroupOrder<T>, public CStorageSet<T>, public COutGroupHandle<T> {
public:
	Generators(uint outGroupMask, char const* name, int degree, int multLenElem = 1) :
		CStorageSet<T>(10, degree * multLenElem),
		COutGroupHandle<T>(outGroupMask, name, degree) {};
	void makeGroupOutput(const CRepository<T>* pElemGroup, bool outToScreen = false, bool checkNestedGroups = true) override;
	virtual void createOrbits(const CRepository<T>* pElemGroup);
	inline auto groupDegree() const			{ return this->m_nc; }
	void resetOrbits() {
		this->setGroupOrder(1);
		this->setStabilizerLengthAut(m_lenStab = groupDegree());
		this->releaseAllObjects();
	}
protected:
	virtual void createOrbitsSet(const CRepository<T>* pElemGroup);
	int testNestedGroups(const CGroupInfo* pElemGroup, CGroupInfo* pRowGroup = NULL, int rowMin = 2, CKOrbits* pKOrb = NULL) const;
	void resetGroup()						{ m_lenStab = groupDegree(); }
	inline bool orbitsCreated() const		{ return m_bOrbitsCreated; }
	inline void setOrbitsCreated(bool val)	{ m_bOrbitsCreated = val; }

private:
	CC virtual bool needUpdate(const T* permRow, const T* pOrbits) {
		if (m_lenStab > this->stabilizerLength())
			return true;
		// Verify if the element is the leading representative of its orbit.
		T i = 0;
		while (permRow[i] == i)
			i++;

		const auto elem = permRow[i];
		return elem == pOrbits[elem];
	}
	CC void savePermutation(ushort degree, const T* permRow, bool rowPermut, bool savePermut) {
		this->addObject(permRow);
		m_lenStab = this->stabilizerLength();
	}

	ushort m_lenStab;
	bool m_bOrbitsCreated = false;
};

void reportNestedGroupCheckResult(int retVal, bool outToScreen);

template<typename T>
void Generators<T>::createOrbits(const CRepository<T>* pElemGroup) {
	// Calculate the orbits and a minimal generating set 
	// of the permutation group under its action on the element set
	if (orbitsCreated())
		return;

	createOrbitsSet(pElemGroup);
	setOrbitsCreated(true);
}

template<typename T>
void Generators<T>::createOrbitsSet(const CRepository<T>* pElemGroup) {
	resetOrbits();

	// Adding orbits:
	auto* pOrb = this->getNextObject();
	for (int i = groupDegree(); i--;)
		pOrb[i] = i;

	// ...  and trivial permutation:
	this->addObject(pOrb);
	const auto groupOrder = pElemGroup->numObjects();
	for (int i = 1; i < groupOrder; i++) {
		const auto* c = pElemGroup->getObject(i);
		this->addAutomorphism(groupDegree(), c, pOrb, true, false, true);
	}

	this->updateGroupOrder(groupDegree(), pOrb);
}

template<typename T>
int Generators<T>::testNestedGroups(const CGroupInfo* pElemGroup, CGroupInfo* pRowGroup, int rowMin, CKOrbits* pKOrb) const {
	const auto pntr = (const alldata*)pElemGroup;
	if (pntr->param(t_nestedGroups) > 1)
		rowMin = 2;
	else
		if (!pRowGroup && !pKOrb)
			return -1;

	tchar ts[MAX_PLAYER_NUMBER];
	ctchar* mi = pntr->result();
	CGroupInfo* pRowGroupOut = NULL;
	const auto groupOrder = pElemGroup->numObjects();
	const auto rowMax = pntr->numDaysResult();
	for (int j = rowMin; j <= rowMax; j++) {
		if (j == rowMax && !(pRowGroupOut = pRowGroup) && !pKOrb)
			break;

		const auto len = pElemGroup->lenObject() * j;
		for (int i = 0; i < groupOrder; i++) {
			pntr->kmSortMatrixForReorderedPlayers(mi, j, pElemGroup->getObject(i), ts, false, pKOrb);
			if (MEMCMP(mi, pntr->transformedMatrix(), len))
				return j;

			if (pRowGroupOut)
				pRowGroupOut->updateGroup(ts);
		}
	}
	return rowMin == rowMax ? -1 : 0;
}

template<typename T>
class RowGenerators : public Generators<T> {
public:
    RowGenerators(uint outGroupMask, int rowNumb, int lenElem = 1) : 
		m_lenElem(lenElem), Generators<T>(outGroupMask, "", rowNumb, lenElem), m_pRowGroup(NULL) {
		m_outMask = 4;
		m_sActionOn = "matrix rows, |Aut(R)|";
		m_bGroupConstructed = false;
	}

	~RowGenerators()						{ delete m_pRowGroup; }
	void makeGroupOutput(const CRepository<T>* pElemInfo, bool outToScreen = false, bool checkNestedGroups = true) override;
	const char* name() override				{ return m_sName.c_str(); }
	virtual int createGroupAndOrbits(const CRepository<tchar>* pElemGroup);
protected:
	int getGroup(const CRepository<tchar>* pElemGroup);
	int lenElem() const						{ return m_lenElem; }
	virtual void createTable(ctchar* pSolution)	{}
	virtual int createGroup(const CRepository<tchar>* pElemGroup) {
		return this->testNestedGroups(pElemGroup, m_pRowGroup, this->lenObject());
	}
	virtual void generatorOutput(bool outToScreen, const char* pErr = NULL) {
		Generators<T>::makeGroupOutput(m_pRowGroup, outToScreen, false);
	}

	CRepository<T>* m_pRowGroup = NULL;
	std::string m_sName;
	std::string m_sActionOn;
	uint m_outMask;
	bool m_bGroupConstructed;
	int m_groupState;
	const int m_lenElem;
};


template<typename T>
void RowGenerators<T>::makeGroupOutput(const CRepository<T>* pElemGroup, bool outToScreen, bool checkNestedGroups) {
#ifndef USE_CUDA
	char errBuf[48], * pErr = NULL;
	const auto retVal = createGroupAndOrbits(pElemGroup);
	if (retVal > 0)
		snprintf(pErr = errBuf, sizeof(errBuf), "Nested groups check failed on row %d\n", retVal);

	if (this->m_outGroupMask & m_outMask) {
		m_sName = std::format("\n{}Orbits and generators of Aut(M) acting on {} = {}",
			(pErr ? pErr : ""), m_sActionOn, m_pRowGroup->orderOfGroup());
		Generators<T>::setName(name());
		generatorOutput(outToScreen, pErr);
		// To avoid writing this error (if any) to the file twice
		pErr = NULL;
	}
	if (this->m_outGroupMask & (m_outMask << 1)) {
		m_sName = std::format("\n{}Aut(M) acting on {} = {}",
			(pErr ? pErr : ""), m_sActionOn, m_pRowGroup->orderOfGroup());
		COutGroupHandle<T>::makeGroupOutput(m_pRowGroup, outToScreen, false);
	}

	m_bGroupConstructed = false;
	reportNestedGroupCheckResult(retVal, outToScreen);
#endif
}


template<typename T>
int RowGenerators<T>::getGroup(const CRepository<tchar>* pElemGroup) {
	if (!m_pRowGroup)
		m_pRowGroup = new CRepository<tchar>(this->lenObject(), 10);
	else
		m_pRowGroup->releaseAllObjects();

	return createGroup(pElemGroup);
}

template<typename T>
int RowGenerators<T>::createGroupAndOrbits(const CRepository<tchar>* pElemGroup) {
	if (m_bGroupConstructed)
		return m_groupState;

	createTable(((alldata*)pElemGroup)->result());
	m_groupState = getGroup(pElemGroup);

	this->createOrbits(m_pRowGroup);
	m_bGroupConstructed = true;
	return m_groupState;
}
#if OUTPUT_VECTOR_STAT
static size_t nMatr = 0;
static size_t nMatrMax = 0;
unsigned char *pMatrixStorage = NULL;
#endif

template<typename T>
void Table<T>::printTable(const T *c, bool outCntr, bool outToScreen, int nl, const char* pStartLine, const int *idx)
{
	char buffer[1024], *pBuf = buffer;
	FOPEN_F(f, this->m_pFileName, "a");
	if (outCntr)
		m_cntr++;

	auto* pName = name();
	if (pName && strlen(pName)) {
		if (outCntr && m_bOutCntr) {
			// When pName starts with '\n', we need to place them before counter.
			while (*pName == '\n') {
				pName++;
				SPRINTFD(pBuf, buffer, "\n");
			}

			SPRINTFD(pBuf, buffer, "%5zd: %s\n", m_cntr, pName);
		}
		else
			SPRINTFD(pBuf, buffer, "%s:\n", pName);

		_printf(f, outToScreen, buffer);
	}

	if (idx) {
		for (int i = 0; i < nl; i++)
			outMatrix(c + idx[i] * m_nc, 1, m_nc, m_np, m_ns, f, m_makeString, outToScreen, pStartLine, -1, NULL, i == nl - 1);
	} else
		outMatrix(c, nl == 0 ? m_nl : nl, m_nc, m_np, m_ns, f, m_makeString, outToScreen, pStartLine);

	outLS((ctchar *)c, nl, f, outToScreen);
#if OUTPUT_VECTOR_STAT
	// Output of a vector with the i-th coordinate equal to the number  
	// of appearances of the i-th player first player in the group.
	if (m_cntr) {
		auto buf = new unsigned char [m_nc];
		memset(buf, 0, m_nc * sizeof(*buf));
		const auto lenMatr = m_nl * nGroups;
		auto matrix = new unsigned char[lenMatr];
		auto pntr = c;
		for (T i = 0; i < m_nl; i++, pntr += m_nc) {
			auto pMatrixRow = matrix + nGroups * i;
			for (T j = 0; j < nGroups; j++)
				buf[pMatrixRow[j] = pntr[j*m_np]]++;
		}

		for (size_t i = 0; i < nMatr; i++) {
			if (memcmp(matrix, pMatrixStorage + lenMatr * i, lenMatr))
				continue;

			delete[] matrix;
			matrix = NULL;
			break;
		}

		if (matrix) {
			// New matrix found
			if (nMatr == nMatrMax) {
				auto pTmp = pMatrixStorage;
				nMatrMax = 2 * (nMatrMax + 1);
				pMatrixStorage = new unsigned char[nMatrMax * lenMatr];
				if (pTmp) {
					memcpy(pMatrixStorage, pTmp, nMatr * lenMatr);
					delete[] pTmp;
				}
			}

			memcpy(pMatrixStorage + nMatr++ * lenMatr, matrix, lenMatr);
			delete[] matrix;
		}

		static char idx[256] = { '\0' };
		pBuf = idx;
		if (idx[0] == '\0') {
			for (T i = 0; i < m_nc; i++)
				SPRINTFD(pBuf, idx, " %2d:", i);

			SPRINTFD(pBuf, idx, "\n");
		}

		pBuf = buffer;
		for (T i = 0; i < m_nc; i++)
			SPRINTFD(pBuf, buffer, " %3d", buf[i]);

		delete[] buf;
		SPRINTFD(pBuf, buffer, "\n");
		_printf(f, false, idx);
		_printf(f, false, buffer);

		if (m_nc == 15 && m_cntr == 101) {
			auto pntr = pMatrixStorage;
			for (size_t i = 0; i < nMatr; i++, pntr += lenMatr) {
				pBuf = buffer;
				SPRINTFD(pBuf, buffer, "\nMatrix #%zd\n", i);
				auto ptr = pntr;
				for (T i1 = 0; i1 < m_nl; i1++, ptr += nGroups) {
					for (T j = 0; j < nGroups; j++)
						SPRINTFD(pBuf, buffer, " %3d", ptr[j]);
					SPRINTFD(pBuf, buffer, "\n");
				}

				_printf(f, false, buffer);
			}

			delete[] pMatrixStorage;
			nMatr = nMatrMax = 0;
		}
	}
#endif
	FCLOSE_F(f);
}

template<typename T>
void Table<T>::printTableInfo(const char* pInfo) {
	if (!pInfo)
		return;

	FOPEN_F(f, this->m_pFileName, "a");
	fprintf(f, pInfo);
	FCLOSE_F(f);
}
