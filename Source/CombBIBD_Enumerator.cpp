#include "CombBIBD_Enumerator.h"
#include "DesignDB.h"

template class CCombBIBD_Enumerator<TDATA_TYPES>;

FClass2(CCombBIBD_Enumerator)::~CCombBIBD_Enumerator() {
	delete[] m_FirstPartSolutionIdx;
	delete[] m_bSolutionsWereConstructed;
	delete m_pSpareMatrix;
	delete[] m_pGroupOrders;
	delete m_pGroupOrder;

	if (m_bColPermutOwner) {
		// Only the master thread is the owner of these data
		delete[] columnPermut();
	}

	delete[] m_pColPermut;
	delete m_pCanonChecker;
}

#if !CONSTR_ON_GPU
FClass2(CCombBIBD_Enumerator, int)::addLambdaInfo(char *buf, size_t lenBuffer, const char* pFormat, size_t *pLambdaSetSize) const {
	return addLambdaInform(paramSet(t_lSet), buf, lenBuffer, pLambdaSetSize);
}

FClass2(CCombBIBD_Enumerator, int)::getJobTitleInfo(char* buffer, int lenBuffer) const {
	const auto inSys = this->getInSys();
	if (inSys->objectType() == t_objectType::t_Kirkman_Triple)
		return SNPRINTF(buffer, lenBuffer, "KirkSys(%" _FRMT", %1"_FRMT ", 1", inSys->rowNumbExt(), inSys->GetK());

	return SNPRINTF(buffer, lenBuffer, "%s(%3" _FRMT", %2" _FRMT", ", getObjName(), inSys->rowNumbExt(), inSys->GetK());
}
#endif

FClass2(CCombBIBD_Enumerator, void)::getEnumerationObjectKey(char* pKey, int len) const {
	char buffer[64];
	this->makeJobTitle(this->designParams(), buffer, countof(buffer));
	SNPRINTF(pKey, len, "%s", buffer + strlen(this->getObjName()));
}

FClass2(CCombBIBD_Enumerator, char*)::getEnumerationObjectKeyA(char* pKey, int len, const char* pKeyIn) {
	if (pKeyIn) {
		if (pKey != pKeyIn)
			strncpy_s(pKey, len, pKeyIn, len - 1);
	}
	else {
		getEnumerationObjectKey(m_pAdjKey = pKey, len);
		m_pAdjKeyPrefix = pKey + len;
	}

	auto* pntr = strstr(pKey, "})");
	*pntr = '\0';

	// Key prefix
	strcpy_s(pKey+len, len, pKey);
	pntr = strstr(pKey+len, "{");
	*pntr = '\0';
	if (!pKeyIn)
		m_lenAdjKey = strlen(pKey);

	return pKey;
}

FClass2(CCombBIBD_Enumerator, int)::compareEnumerationDB_record(const char* record) {
	char keyCmp[64];
	const auto len = countof(keyCmp) / 2;
	char* keyCmpPrefix = keyCmp + len;
	getEnumerationObjectKeyA(keyCmp, len, record);
	auto resCmp = strcmp(keyCmpPrefix, m_pAdjKeyPrefix);
	if (!resCmp) {
		const auto lenKeyCmp = strlen(keyCmp);
		if (lenKeyCmp < m_lenAdjKey)
			resCmp = -1;
		else
			resCmp = lenKeyCmp == m_lenAdjKey ? strcmp(keyCmp, m_pAdjKey) : 1;
	}
	return resCmp;
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
	// Create matrix which will be used only when there is non-trivial group on parts
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

FClass2(CCombBIBD_Enumerator, void)::CreateAuxiliaryStructures(EnumeratorPntr pMaster) {
	if (m_pCanonChecker ||							// canonicity checker was already constructed
		!designParams()->find_master_design ||      // no need to find master design OR
		!pMaster && designParams()->create_commonData()) // using threads, but now we are in master   
		return;

	const auto b = matrix()->colNumb();
	const auto v = matrix()->rowNumb() - 1;
	auto pCombBIBD_master = static_cast<const CCombBIBD_Enumerator*>(pMaster);
	setMaster(pMaster);
	if (pMaster) {
		m_pColumnPermut = (T*)pCombBIBD_master->columnPermut();
		if (m_pColumnPermut)
			m_pColPermut = new T[b];

		initDesignDB(pMaster, 1);
	}

	// Creating structures for the search of "master" designs for Combined BIBDs
	auto *pOriginalMatrix = new C_InSys<T, S>();
	pOriginalMatrix->Init(v, b);

	m_pCanonChecker = new CMatrixCanonChecker<TDATA_TYPES>(pOriginalMatrix, t_matrixOwner);
	m_pCanonChecker->initiateColOrbits(v, 0, NULL, true);

	T lambda = 0;
	if (getInSys()->objectType() != t_objectType::t_Kirkman_Triple) {
		// Set first two rows and first columns of "master" design
		// Because they always be the same, it would be better to do it only once
		// First, we need to define parameter lambda for "master" design
		const auto lambdaSet = paramSet(t_lSet);

		for (auto j = numParts(); j--;)
			lambda += lambdaSet->GetAt(j);
	}
	else {
		lambda = 1;
		m_nNum_lambdas = 2; // their values are 0 and 1
	}

	const auto k = this->getInSys()->GetNumSet(t_kSet)->GetAt(0);
	const auto r = static_cast<T>(lambda * (v - 1) / (k - 1));
	T solution[] = { r, lambda, static_cast<T>(r - lambda) };
	m_pCanonChecker->MakeRow(0, solution);
	m_pCanonChecker->MakeRow(1, solution+1);
	auto* pRow = pOriginalMatrix->GetRow(1);
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

	const auto v = matrix()->rowNumb() - 1;
	const auto k = getInSys()->GetNumSet(t_kSet)->GetAt(0);
	const auto pR_set = paramSet(t_rSet);
	const auto lambdaSet = paramSet(t_lSet);

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

	if (designParams()->find_master_design)
		initDesignDB(NULL, 1);   // Create DB for storing "master" BIBDs
}

FClass2(CCombBIBD_Enumerator, void)::ConstructedDesignProcessing() const {
	// Merging parts into one BIBD with the first two canonical rows.
	const auto v = matrix()->rowNumb() - 1;
	const auto b = matrix()->colNumb();
	// Creating:
	// (a) "original" matrix by merging parts according to columnPermut()
	// (b) the orbits on columns for each row of that matrix

	// Column permutation can be changed later, so copy it.
	T colPermut[256];
	T* pColPermut = (countof(colPermut) >= b) ? colPermut : new T[b];
	memcpy(pColPermut, columnPermut(), b * sizeof(pColPermut[0]));
	auto* pMatr = m_pCanonChecker->matrix();
	for (T i = 2; i < v; i++) {
		auto* pRowSrc = matrix()->GetRow(i + 1);
		auto* pRow = pMatr->GetRow(i);
		for (auto j = b; --j >= 1;)
			pRow[j] = pRowSrc[pColPermut[j]];

		m_pCanonChecker->CreateColumnOrbits(i, pRow, pColPermut);
	}

	if (pColPermut != colPermut)
		delete[] pColPermut;

	m_pCanonChecker->CanonizeMatrix();
	AddMatrixToDB(m_pCanonChecker, 1);
}

FClass2(CCombBIBD_Enumerator, void)::beforeEnumInfoOutput() const {
	if (!designDB())
		return;

	// Sorting "master" BIBDs by their numbers of decompositions
	outBlockTitle("Decomposition info", false);
	designDB()->SortRecods(outFile());
	outString(" \n" END_OUT_BLOCK "Decomposition info " BEG_OUT_BLOCK "\n\n", outFile());
}