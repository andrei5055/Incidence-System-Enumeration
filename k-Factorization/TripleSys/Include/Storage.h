#pragma once

template<typename T>
class CStorage {
public:
	CC CStorage(int numObjects, int lenObj = 1) : m_lenObj(lenObj) {
		if (numObjects)
			m_pObjects = new T[numObjects * m_lenObj];
	}
	CC ~CStorage()								{ delete[] m_pObjects; }
	CC T* reallocStorageMemory(int numObjects)	{
		return ::reallocStorageMemory(&m_pObjects, m_lenObj * numObjects);
	}
	CC virtual T* getObject(int idx = 0) const	{ return m_pObjects + idx * m_lenObj; }
	CC inline T** getObjectsPntr()				{ return &m_pObjects; }
	CC inline auto lenObject() const			{ return m_lenObj; }
	CC inline void setLenCompare(int len)		{ m_lenCompare = len; }
	CC T& operator[](int i) const				{ return *getObject(i); }
	CC uint findObject(const T *pObj, uint low, uint high) const {
		// search for element
		high--;
		while (low <= high) {
			const auto itr = low + ((high - low) >> 1);
			const auto cmp = compareObjects(itr, pObj);
			if (!cmp)
				return itr;

			if (low == high)
				break;

			if (cmp < 0)
				low = itr + 1;  // ignore left half
			else
				high = itr - 1; // ignore right half
		}

		return UINT_MAX;		// not found 
	}
	CC int getElementIndex(const T* tr, uint nElem) const {
		// search for element 	
		int itr;
		int low = 0;
		auto high = itr = nElem - 1;
		int cmp = -1;
		while (low <= high) {
			itr = low + ((high - low) >> 1);
			cmp = MEMCMP(getObject(itr), tr, m_lenObj);
			if (!cmp)
				return -itr - 1;

			if (cmp < 0)
				low = itr + 1;  // ignore left half
			else
				high = itr - 1; // ignore right half
		}

		if (cmp < 0)
			itr++;

		return itr;
	}
protected:
	const int m_lenObj;
private:
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
	return (int)(*getObject(idx)) - *obj;
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
		::reallocStorageMemory(&m_pIdx, CStorageSet<T>::m_numObjectsMax * sizeof(m_pIdx[0]));
		return retVal;
	}
	CC void push_back(int idx) { m_pIdx[CStorageSet<T>::m_numObjects++] = idx; }
	CC void insert(int idx, int value) {
		const auto len = (CStorageSet<T>::m_numObjects++ - idx) * sizeof(m_pIdx[0]);
		const auto pSrc = m_pIdx + idx;
		MEMMOVE((void*)(pSrc + 1), pSrc, len);
		m_pIdx[idx] = value;
	}
	CC T* getObjAddr(int idx) { return idx < CStorageSet<T>::m_numObjectsMax ? CStorage<T>::getObject(idx) : reallocStorageMemory(); }
	CC inline int* getIndices() const { return m_pIdx; }
	CC CStorageIdx& operator = (const CStorageIdx& other) {
		auto** ppStorage = CStorage<T>::getObjectsPntr();
		delete[] * ppStorage;
		delete[] m_pIdx;
		m_pIdx = NULL;	// We don't need indices for now 	
		CStorageSet<T>::m_numObjectsMax = CStorageSet<T>::m_numObjects = other.numObjects();
		const auto len = CStorageSet<T>::m_numObjects * CStorage<T>::m_lenObj;
		*ppStorage = new T[len];
		memcpy(*ppStorage, other.CStorage<T>::getObject(), len);
		return *this;
	}
	CC void copyIndex(const CStorageIdx& other) {
		m_pIdx = new int[other.numObjects()];
		memcpy(m_pIdx, other.getIndices(), CStorageSet<T>::m_numObjects * sizeof(m_pIdx[0]));
	}
/*
	CC int getElementIndex(const T* tr, uint nElem) const {
		// search for element 	
		int itr;
		int low = 0;
		auto high = itr = nElem - 1;
		int cmp = -1;
		while (low <= high) {
			itr = low + ((high - low) >> 1);
			cmp = MEMCMP(getObject(itr), tr, CStorage<T>::lenObject());
			if (!cmp)
				return -itr - 1;

			if (cmp < 0)
				low = itr + 1;  // ignore left half
			else
				high = itr - 1; // ignore right half
		}

		if (cmp < 0)
			itr++;

		return itr;
	}
	*/
private:
	int* m_pIdx = NULL;
};

