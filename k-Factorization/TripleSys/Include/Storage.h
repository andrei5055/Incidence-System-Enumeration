#pragma once

template<typename T>
class CStorage {
public:
	CC CStorage(int numObjects, int lenObj = 1) : m_lenObj(lenObj) {
		setLenCompare(lenObj);
		if (numObjects)
			m_pObjects = new T[numObjects * m_lenObj];
	}
	CC ~CStorage()								{ delete[] m_pObjects; }
	CC inline T* reallocStorageMemory(int numObjects) {
		return ::reallocStorageMemory(&m_pObjects, m_lenObj * numObjects);
	}
	CC virtual T* getObject(int idx = 0) const	{ return m_pObjects + idx * m_lenObj; }
	CC inline T** getObjectsPntr()				{ return &m_pObjects; }
	CC inline auto lenObject() const			{ return m_lenObj; }
	CC inline void setLenCompare(int len)		{ m_lenCompare = len; }
	CC inline T& operator[](int i) const		{ return *getObject(i); }
	CC inline uint findObject(const T* pObj, uint low, uint high) const {
		return static_cast<uint>(binarySearch(pObj, low, high, false));
	}
	CC inline int getElementIndex(const T* obj, uint nElem) const {
		return binarySearch(obj, 0, nElem, true);
	}
protected:
	const int m_lenObj;
private:
	CC int binarySearch(const T* pObj, long low, long high, bool returnInsertionPoint) const {
		if (high == 0 || low > high)  // Guard against invalid ranges
			return returnInsertionPoint && !high ? 0 : UINT_MAX;

		high--;
		long mid = 0;
		int cmp = -1;
		while (low <= high) {
			mid = low + ((high - low) >> 1);
			cmp = compareObjects(mid, pObj);

			if (cmp == 0)
				return returnInsertionPoint ? -mid - 1 : mid;

			if (cmp < 0)
				low = mid + 1;
			else
				high = mid - 1;
		}

		// Not found
		if (returnInsertionPoint)
			return (cmp < 0) ? mid + 1 : mid;
		else
			return UINT_MAX;
	}

	CC int compareObjects(uint idx, const T* obj) const;
	T* m_pObjects = NULL;
	int m_lenCompare;
};

template<>
inline int CStorage<tchar>::compareObjects(uint idx, const tchar* obj) const {
	return MEMCMP(getObject(idx), obj, m_lenCompare);
}

template<>
inline int CStorage<uint>::compareObjects(uint idx, const uint *obj) const {
	const auto val = *getObject(idx);
	return val > *obj? 1 : (val < *obj? -1 : 0);
}

template<typename T>
class CStorageSet : public CStorage<T> {
public:
	CC CStorageSet(int numObjects, int lenObj = 1) : CStorage<T>(numObjects, lenObj) {
		m_numObjectsMax = numObjects;
	}
	CC inline auto numObjects() const		{ return m_numObjects; }
	CC inline auto numObjectsMax() const	{ return m_numObjectsMax; }
	CC inline auto getNextObject()          { 
		if (m_numObjects >= m_numObjectsMax)
			reallocStorageMemory();
		return CStorage<T>::getObject(m_numObjects++);
	}

	CC inline auto* addObject(const T* pObj) {
		return (const T*)memcpy(getNextObject(), pObj, this->lenObject());
	}
	CC inline void releaseObject()			{ m_numObjects--; }
	CC inline void releaseAllObjects()		{ m_numObjects = 0; }
	int m_numObjects = 0;
protected:
	CC T* reallocStorageMemory() {
		return CStorage<T>::reallocStorageMemory(m_numObjectsMax <<= 1);
	}

	int m_numObjectsMax = 0;
};

template<typename T>
class CStorageIdx : public CStorageSet<T> {
public:
	CC CStorageIdx(int numObjects, int lenObj = 1) : CStorageSet<T>(numObjects, lenObj) {
		if (numObjects)
			m_pIdx = new int[numObjects];
	}
	CC ~CStorageIdx() { delete[] m_pIdx; }
	CC inline T* getObjPntr(int idx) const { return CStorage<T>::getObject(m_pIdx[idx]); }
	CC virtual T* getObject(int idx = 0) const { return CStorage<T>::getObject(m_pIdx[idx]); }
	CC T* reallocStorageMemory() {
		const auto retVal = CStorageSet<T>::reallocStorageMemory();
		::reallocStorageMemory(&m_pIdx, this->m_numObjectsMax * sizeof(m_pIdx[0]));
		return retVal;
	}
	CC inline void push_back(int idx) { m_pIdx[this->m_numObjects++] = idx; }
	CC inline void insert(int idx, int value) {
		const auto pSrc = m_pIdx + idx;
		MEMMOVE((void*)(pSrc + 1), pSrc, (this->m_numObjects++ - idx) * sizeof(m_pIdx[0]));
		m_pIdx[idx] = value;
	}
	CC T* getObjAddr(int idx) { return idx < this->m_numObjectsMax ? CStorage<T>::getObject(idx) : reallocStorageMemory(); }
	CC inline int* getIndices() const { return m_pIdx; }
	CC CStorageIdx& operator = (const CStorageIdx& other) {
		auto** ppStorage = CStorage<T>::getObjectsPntr();
		delete[] * ppStorage;
		delete[] m_pIdx;
		m_pIdx = NULL;	// We don't need indices for now 	
		this->m_numObjectsMax = this->m_numObjects = other.numObjects();
		const auto len = this->numObjects() * CStorage<T>::m_lenObj;
		*ppStorage = new T[len];
		memcpy(*ppStorage, other.CStorage<T>::getObject(), len);
		return *this;
	}
	CC inline void copyIndex(const CStorageIdx& other) {
		m_pIdx = new int[other.numObjects()];
		memcpy(m_pIdx, other.getIndices(), this->numObjects() * sizeof(m_pIdx[0]));
	}
	CC inline auto isProcessed(ctchar* tr) {
		const auto numRegisteredTrs = this->numObjects();
		updateRepo(tr);
		return numRegisteredTrs == this->numObjects();
	}
	CC int updateRepo(const T* tr) {
		// search for element 
		const auto nElem = this->numObjects();
		const auto itr = CStorageSet<T>::getElementIndex(tr, nElem);
		if (itr < 0)
			return itr;

		auto* cmpTr = getObjAddr(nElem);

		if (itr < nElem)
			insert(itr, nElem);
		else
			push_back(nElem);

		memcpy(cmpTr, tr, this->lenObject());
		return itr;
	}

private:
	int* m_pIdx = NULL;
};

