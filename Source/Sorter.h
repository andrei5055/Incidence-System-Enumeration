#pragma once


template<typename T>
class CSorter {
	typedef int (*compareFunc)(const T *pRec1, const T *pRec2);
	typedef int (*compareFuncA)(const T* pRec1, const T* pRec2, const CSorter *pntr);
public:
	CSorter(size_t len)								{ setRecordLength(len); }
	~CSorter()										{ setRecordStorage(); }
	inline size_t recordLength() const				{ return m_nRecLen; }
protected:
	void quickSort(size_t* arr, size_t left, size_t right) const;
	const T *getRecord(size_t idx) const			{ return firstRecord() + idx * recordLength(); }
	inline const T *firstRecord() const				{ return m_pRecStorage; }
	inline void setRecordLength(size_t len)			{ m_nRecLen = len; }
	inline void setRecordStorage(const T *pntr = NULL)	{ if (m_bStorageOwner) delete[] m_pRecStorage;
													  m_pRecStorage = pntr;
													}
	inline void setCompareFuncA(compareFuncA pF)	{ m_pCompareFuncA = pF; }
	inline void setCompareFunc(compareFunc pF)		{ m_pCompareFunc = pF; }
	inline void setStorageOwner(bool val = true)	{ m_bStorageOwner = val; }
private:
	size_t m_nRecLen;							// length of each record
	const T *m_pRecStorage = NULL;					// memory for record storing
	compareFunc m_pCompareFunc = NULL;
	compareFuncA m_pCompareFuncA = NULL;
	bool m_bStorageOwner = false;
};

template<typename T>
void CSorter<T>::quickSort(size_t* arr, size_t left, size_t right) const {
	size_t i = left, j = right;
	const auto pivotIdx = (left + right) >> 1;
	auto pivot = getRecord(arr[pivotIdx]);

	/* partition */
	while (i <= j) {
		if (m_pCompareFuncA) {
			while (i != pivotIdx && m_pCompareFuncA(getRecord(arr[i]), pivot, this) == -1)
				i++;

			while (j != pivotIdx && m_pCompareFuncA(getRecord(arr[j]), pivot, this) == 1)
				j--;
		} else{
			while (i != pivotIdx && m_pCompareFunc(getRecord(arr[i]), pivot) == -1)
			i++;

			while (j != pivotIdx && m_pCompareFunc(getRecord(arr[j]), pivot) == 1)
				j--;
		}

		if (i < j) {
			const auto tmp = arr[i];
			arr[i++] = arr[j];
			arr[j--] = tmp;
		}
		else {
			if (i == j)
				i++;
			break;
		}
	}

	/* recursion */
	if (left < j)
		quickSort(arr, left, j);

	if (i < right)
		quickSort(arr, i, right);
}