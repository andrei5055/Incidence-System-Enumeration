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
	CC inline T* getObject(int idx = 0) const	{ return m_pObjects + idx * m_lenObj; }
	CC inline T** getObjectsPntr()				{ return &m_pObjects; }
	CC inline auto lenObject() const			{ return m_lenObj; }
protected:
	const int m_lenObj;
private:
	T* m_pObjects = NULL;
};

template<typename T>
class CStorageIdx : public CStorage<T> {
public:
	CC CStorageIdx(int numObjects, int lenObj = 1) : CStorage<T>(numObjects, lenObj) {
		if (numObjects)
			m_pIdx = new int[m_numObjectsMax = numObjects];
	}
	CC ~CStorageIdx() { delete[] m_pIdx; }
	CC inline T* getObjPntr(int idx) const { return CStorage<T>::getObject(m_pIdx[idx]); }
	CC inline auto numObjects() const { return m_numObjects; }
	CC T* reallocStorageMemory() {
		m_numObjectsMax <<= 1;
		::reallocStorageMemory(&m_pIdx, m_numObjectsMax * sizeof(m_pIdx[0]));
		return CStorage<T>::reallocStorageMemory(m_numObjectsMax);
	}
	CC void push_back(int idx) { m_pIdx[m_numObjects++] = idx; }
	CC void insert(int idx, int value) {
		const auto len = (m_numObjects++ - idx) * sizeof(m_pIdx[0]);
		const auto pSrc = m_pIdx + idx;
		MEMMOVE((void*)(pSrc + 1), pSrc, len);
		m_pIdx[idx] = value;
	}
	CC T* getObjAddr(int idx) { return idx < m_numObjectsMax ? CStorage<T>::getObject(idx) : reallocStorageMemory(); }
	CC inline int* getIndices() const { return m_pIdx; }
	CC CStorageIdx& operator = (const CStorageIdx& other) {
		auto** ppStorage = CStorage<T>::getObjectsPntr();
		delete[] * ppStorage;
		delete[] m_pIdx;
		m_pIdx = NULL;	// We don't need indices for now 	
		m_numObjectsMax = m_numObjects = other.numObjects();
		const auto len = m_numObjects * CStorage<T>::m_lenObj;
		*ppStorage = new T[len];
		memcpy(*ppStorage, other.getObject(), len);
		return *this;
	}
	CC void copyIndex(const CStorageIdx& other) {
		m_pIdx = new int[other.numObjects()];
		memcpy(m_pIdx, other.getIndices(), m_numObjects * sizeof(m_pIdx[0]));
	}

	int m_numObjects = 0;
protected:
	int m_numObjectsMax = 0;
private:
	int* m_pIdx = NULL;
};

