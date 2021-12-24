#include "CombBIBD_Enumerator.h"

template class CCombBIBD_Enumerator<TDATA_TYPES>;

#if !CONSTR_ON_GPU
FClass2(CCombBIBD_Enumerator, int)::addLambdaInfo(char *buf, size_t lenBuffer, const char* pFormat, size_t *pLambdaSetSize) const {
	return addLambdaInform(paramSet(t_lSet), buf, lenBuffer, pLambdaSetSize);
}

FClass2(CCombBIBD_Enumerator, int)::getJobTitleInfo(char* buffer, int lenBuffer) const {
	const auto inSys = this->getInSys();
	return SNPRINTF(buffer, lenBuffer, "%s(%3" _FRMT", %2" _FRMT", ", getObjName(), inSys->rowNumbExt(), inSys->GetK());
}
#endif

FClass2(CCombBIBD_Enumerator, void)::getEnumerationObjectKey(char* pKey, int len) const {
	char buffer[64];
	this->makeJobTitle(this->designParams(), buffer, countof(buffer));
	SNPRINTF(pKey, len, "%s", buffer + strlen(this->getObjName()));
}

FClass2(CCombBIBD_Enumerator, char*)::getEnumerationObjectKeyA(char* pKey, int len) const {
	getEnumerationObjectKey(pKey, len);
	auto* pntr = strstr(pKey, "})");
	*pntr = '\0';
	return pKey;
}

FClass2(CCombBIBD_Enumerator, RowSolutionPntr)::setFirstRowSolutions() {
	auto pSolutions = this->rowStuff(1);
	const auto* pR_set = this->paramSet(t_rSet);
	for (auto i = this->numParts(); i--;) {
		auto pPartSolutions = pSolutions + i;
		pPartSolutions->InitSolutions(1, 1);
		pPartSolutions->AddElement(pR_set->GetAt(i));
	}

	return pSolutions;
}

FClass2(CCombBIBD_Enumerator, void)::CreateForcedRows() {
	// As of today (12/12/2021), method is called for master only
	// Combined BIBDs having n parts will be constructed with first "artificial" row:
	//  n, n, ...,n, n-1, ..., n-1, n-2,..., n-2,... 2,...2,1...,1
	// It's why we need to start enumeration with the row number 1.
	this->setCurrentRowNumb(1);
	CreateFirstRow();
	if (designParams()->find_master_design)
		createColumnPermut();
}

FClass2(CCombBIBD_Enumerator, void)::CreateFirstRow() {
	const auto pInSys = this->getInSys();
	const auto k = pInSys->GetNumSet(t_kSet)->GetAt(0);
	const auto v = pInSys->rowNumb() - 1;
	const auto pR_set = this->paramSet(t_rSet);
	const auto iMax = static_cast<T>(this->numParts());
	S *pFirstRow = pInSys->GetRow(0);

	for (T i = 0; i < iMax; i++) {
		const auto b = pR_set->GetAt(i) * v / k;
		rowSetFragm(pFirstRow, iMax - i, b);
		pFirstRow += b;
	}
}

FClass2(CCombBIBD_Enumerator, CMatrixData<T, S> *)::CreateSpareMatrix(const EnumeratorPntr pMaster) {
	// Create matris which will be used only when there is non-trivial group on parts
	if (m_pSpareMatrix)
		return m_pSpareMatrix;

	if (!pMaster)
		pMaster = this;

	const auto pMatr = pMaster->matrix();
	MatrixDataPntr pSpareMatrix = m_pSpareMatrix = new CMatrixData<T, S>();
	pSpareMatrix->Init(pMatr->rowNumb(), pMatr->colNumb());
	memcpy(pSpareMatrix->GetRow(0), pMatr->GetRow(0), pMatr->colNumb() * sizeof(pSpareMatrix->GetRow(0)[0]));
	auto *pPartsInformation = pSpareMatrix->InitPartsInfo(this->numParts());
	pPartsInformation->CopyPartInfo(pMatr->partsInfo());
	return pSpareMatrix;
}

FClass2(CCombBIBD_Enumerator, CGroupOnParts<T> *)::makeGroupOnParts(const CanonicityCheckerPntr owner) {
	const auto lambdaSet = paramSet(t_lSet);
	const auto jMax = lambdaSet->GetSize() - 1;
	CVector<T> lengths;
	T prevLambda = 0;
	uint count;
	uint factorial;
	for (int j = 0; j <= jMax; j++) {
		const auto lambda = lambdaSet->GetAt(j);
		if (prevLambda == lambda) {
			factorial *= (++count);
			if (j < jMax)
				continue;
			j++;
		}

		if (prevLambda) {
			if (factorial > 1) {
				lengths.AddElement(j - count); 	// index of the first BIBD with the same lambda
				lengths.AddElement(count);		// number of BIBDs with the same lambda
				lengths.AddElement(factorial);  // group order
			}
		}
		prevLambda = lambda;
		factorial = count = 1;
	}

	if (!lengths.GetSize())
		return NULL;


	updateCanonicityChecker(rowNumb(), colNumb());
	auto pGroupOnParts = new CGroupOnParts<T>(owner, lengths, 3);
	InitGroupOderStorage(pGroupOnParts);
	return pGroupOnParts;
}

FClass2(CCombBIBD_Enumerator, void)::CreateAuxiliaryStructures(const EnumeratorPntr pMaster) {
	if (m_pCanonChecker ||							// canonicity checker was already constructed
		!designParams()->find_master_design ||      // no need to find master design OR
		designParams()->threadNumb && !pMaster)     // using threads, but now we are in master
		return;


	const auto b = matrix()->colNumb();
	if (pMaster) {
		m_pColumnPermut = (T*)(static_cast<const CCombBIBD_Enumerator*>(pMaster)->columnPermut());
		if (m_pColumnPermut)
			m_pColPermut = new T[b];
	}

	// Creating structures for the search of "master" designs for Combined BIBDs
	const auto v = matrix()->rowNumb() - 1;
	auto *pOriginalMatrix = new CMatrixData<T, S>();
	pOriginalMatrix->Init(v, b);

	m_pCanonChecker = new CMatrixCanonChecker<TDATA_TYPES>(pOriginalMatrix, t_matrixOwner);
	m_pCanonChecker->initiateColOrbits(v, 0, NULL, /*matrix()->partsInfo(),*/ false);// , use_master_sol, pMaster);

	// Set first two rows and first columns of "master" design
	// Because they always be the same, it would be better to do it only once
	// First, we need to define paramete lambda for "master" design
	const auto lambdaSet = this->paramSet(t_lSet);
	T lambda = 0;
	for (auto j = numParts(); j--;)
		lambda += lambdaSet->GetAt(j);

	const auto pInSys = this->getInSys();
	const auto k = pInSys->GetNumSet(t_kSet)->GetAt(0);
	const auto r = static_cast<T>(lambda * (v - 1) / (k - 1));
	T solution[] = { r, lambda, static_cast<T>(r - lambda) };
	m_pCanonChecker->MakeRow(0, solution);
	m_pCanonChecker->MakeRow(1, solution+1);
	return;
	auto* pRow = pOriginalMatrix->GetRow(0);
	rowSetFragm(pRow, 1, r);				// 1's of the first r column of the first row
	rowSetFragm(pRow + r, 0, 2*b - r);		// 0's to remaining part of the first and the whole second row
	rowSetFragm(pRow += b, 1, lambda);		// 1's to first lambda columns of the second row
	rowSetFragm(pRow + r, 1, r - lambda);   // 1's to columns [r, 2*r - lambda - 1] of the second row

	// Remaining part of the first column
	T i = 1;
	while (++i < k)
		*(pRow += b) = 1;

	while (i++ < v)
		*(pRow += b) = 0;
}

FClass2(CCombBIBD_Enumerator, void)::createColumnPermut() {
	// Construct permutation of columns used by all threads when find_master_design is set to 1
	const auto b = matrix()->colNumb();
	m_pColumnPermut = new T[b];
	m_bColPermutOwner = true;		 // only master thread could be the owner
	if (!designParams()->threadNumb) // if the threads will not be used
		m_pColPermut = new T[b];

	const auto pInSys = this->getInSys();
	const auto k = pInSys->GetNumSet(t_kSet)->GetAt(0);
	const auto v = pInSys->rowNumb() - 1;
	const auto pR_set = this->paramSet(t_rSet);
	const auto lambdaSet = this->paramSet(t_lSet);

	T nextIdx, idx, n, i = 0;
	// Indices of first lambda columns of all parts
	// they have 1's in first & second rows
	for (T j = nextIdx = 0; j < numParts(); j++) {
		const auto nMax = lambdaSet->GetAt(j);
		idx = nextIdx + (n = 0);
		while (n++ < nMax)
			m_pColumnPermut[i++] = idx++;

		nextIdx += pR_set->GetAt(j) * v / k;
	}

	// Indices of columns which have 1's in first row, but not in second one
	for (T j = nextIdx = 0; j < numParts(); j++) {
		const auto nMax = pR_set->GetAt(j);
		idx = nextIdx + (n = lambdaSet->GetAt(j));
		while (n++ < nMax)
			m_pColumnPermut[i++] = idx++;

		nextIdx += nMax * v / k;
	}

	// Indices of columns which have 1's in second row, but not in first one
	for (T j = nextIdx = 0; j < numParts(); j++) {
		const auto nMax = 2 * (n = pR_set->GetAt(j)) - lambdaSet->GetAt(j);
		idx = nextIdx + n;
		while (n++ < nMax)
			m_pColumnPermut[i++] = idx++;

		nextIdx += pR_set->GetAt(j) * v / k;
	}

	// Indices of columns which have 0's in first & second rows
	for (T j = nextIdx = 0; j < numParts(); j++) {
		const auto r = pR_set->GetAt(j);
		const auto nMax = r * v / k;
		idx = nextIdx + (n = 2 * r - lambdaSet->GetAt(j));
		while (n++ < nMax)
			m_pColumnPermut[i++] = idx++;

		nextIdx += nMax;
	}
}


FClass2(CCombBIBD_Enumerator, void)::FindMasterBIBD() {
#define USE_COMBINED_MATRIX   false // (!TEST || 1)  // Cannonicity is checked on combined matrix,
											// which has indices of parts in its first row

	// Merging parts into one BIBD with the first two canonical rows.
	const auto v = matrix()->rowNumb() - 1;
	const auto b = matrix()->colNumb();
	// Creating:
	// (a) "original" matrix by merging parts according to columnPermut()
	// (b) the orbits on columns for each row of that matrix
	auto* pMatr = m_pCanonChecker->matrix();
	for (T i = 2; i < v; i++) {
		auto* pRowSrc = matrix()->GetRow(i + 1);
		auto* pRow = pMatr->GetRow(i);
		for (T j = 1; j < b; j++)
			pRow[j] = pRowSrc[columnPermut()[j]];
	}

#if TEST
	pMatr->printOut(this->outFile(), v, 0, this);
#endif

	// Createin orbits on colums

	T nPart = 1;
	T level;
#if USE_COMBINED_MATRIX
	memcpy(m_pColPermut, columnPermut(), b * sizeof(m_pColPermut[0]));
#else
	for (auto j = b; j--;)
		m_pColPermut[j] = j;
#endif
	TestCanonParams<T, S> canonParam =
#if USE_COMBINED_MATRIX
			{ this, matrix(), 1, &nPart, &level, NULL, NULL, m_pColPermut, 1 };
#else
	        { this, pMatr, 1, &nPart, &level, NULL, NULL, m_pColPermut};
#endif
	auto canonMatrix = this->TestCanonicity(v + USE_COMBINED_MATRIX, &canonParam);
}
