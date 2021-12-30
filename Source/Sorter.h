#pragma once


template<typename T>
class CSorter {
	typedef int (*compareFunc)(const T *pRec1, const T *pRec2);
public:
	CSorter(size_t len)							{ setRecordLength(len); }
	~CSorter()									{ delete[] m_pRecStorage; }
protected:
	void quickSort(size_t* arr, size_t left, size_t right) const;
	T *getRecord(size_t idx) const				{ return firstRecord() + idx * recordLength(); }
	inline T *firstRecord() const				{ return m_pRecStorage; }
	inline void setRecordLength(size_t len)		{ m_nRecLen = len; }
	inline size_t recordLength() const			{ return m_nRecLen; }
	inline void setRecordStorage(T *pntr)		{ delete[] m_pRecStorage; m_pRecStorage = pntr; }
	inline void setCompareFunc(compareFunc pFunc) { m_pCompareFunc = pFunc; }
private:
	size_t m_nRecLen;							// length of each record
	T *m_pRecStorage = NULL;					// memory for record storing
	compareFunc m_pCompareFunc;
};

template<typename T>
void CSorter<T>::quickSort(size_t* arr, size_t left, size_t right) const {
	size_t i = left, j = right;
	const auto pivotIdx = (left + right) >> 1;
	auto pivot = firstRecord() + arr[pivotIdx] * recordLength();

	/* partition */
	while (i <= j) {
		while (i != pivotIdx && m_pCompareFunc(getRecord(arr[i]), pivot) == -1)
			i++;

		while (j != pivotIdx && m_pCompareFunc(getRecord(arr[j]), pivot) == 1)
			j--;

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