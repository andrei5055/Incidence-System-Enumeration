#include "BIBD_Enumerator.h"

template class CBIBD_Enumerator<MATRIX_ELEMENT_TYPE>;

template<class T>
int CBIBD_Enumerator<T>::unforcedElement(const CColOrbit<T> *pOrb, int nRow) const
{
	const size_t diffWeight = this->getInSys()->GetK() - pOrb->columnWeight();
	if (diffWeight)
		return diffWeight == this->rowNumb() - nRow ? 1 : -1;

	// All units are there
	return 0;
}

#if !CONSTR_ON_GPU
template<class T>
bool CBIBD_Enumerator<T>::makeFileName(char *buffer, size_t lenBuffer, const char *ext) const
{
	const auto dirLength = this->getDirectory(buffer, lenBuffer) ;
	SNPRINTF(buffer + dirLength, lenBuffer - dirLength, ME_FRMT"_" ME_FRMT"_" ME_FRMT"%s", this->rowNumb(),
						this->getInSys()->GetK(), this->getInSys()->lambda(), ext ? ext : FILE_NAME(""));
	return true;
}

template<class T>
bool CBIBD_Enumerator<T>::makeJobTitle(char *buffer, int lenBuffer, const char *comment) const
{
	const auto v = this->rowNumb();
	const auto b = this->matrix()->colNumb();
	const auto k = this->getInSys()->GetK();
	auto len = SNPRINTF(buffer, lenBuffer, "%s(%3" _FRMT", %3" _FRMT", %2" _FRMT", %2" _FRMT", ", getObjName(), v, b, b * k / v, k);
	len += addLambdaInfo(buffer + len, lenBuffer - len);
	SNPRINTF(buffer + len, lenBuffer - len, ")%s", comment);
	return true;
}
#endif
