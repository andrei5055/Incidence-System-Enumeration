//
//  CanonicityChecker.cpp
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 3/28/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#include "Enumerator.h"

template class CCanonicityChecker<TDATA_TYPES>;

CanonicityChecker(void)::InitCanonicityChecker(T nRow, T nCol, int rank, char *pMem)
{
	m_pPermutRow = (CPermut *)(pMem += nRow);
	m_pPermutRow->Init(nRow, (T *)(pMem += sizeof(*m_pPermutRow)));
	m_pPermutCol = (CPermut *)(pMem += nRow);
	m_pPermutCol->Init(nCol, (T *)(pMem += sizeof(*m_pPermutCol)));
	setColIndex((T *)(pMem += nCol * sizeof(T)));
	m_pCounter = (CCounter<int> *)(pMem += (2 * nCol * sizeof(T)));
	m_pCounter->Init(rank, (int *)(pMem += sizeof(*m_pCounter)));
	setPermStorage((PermutStoragePntr)(pMem += rank * sizeof(int)));
	m_nColNumbStorage = (CColNumbStorage **)(pMem += sizeof(Class2(CPermutStorage)));
	pMem += rank * sizeof(*m_nColNumbStorage);
	for (int i = rank; i--;) {
		m_nColNumbStorage[i] = (CColNumbStorage *)(pMem);
		m_nColNumbStorage[i]->Init(nCol, (T *)(pMem += sizeof(*m_nColNumbStorage[i])));
		pMem += nCol * sizeof(T);
	}

	m_pImprovedSol = NULL;
	setSolutionStorage(NULL);
}

CanonicityChecker(T *)::init(T nRow, T numParts, bool savePerm, T *pOrbits, T** pPermRows, bool groupOnParts, T *pPermCol) {
	T *pRow, *pCol;
	if (!pPermCol) {
		if (!groupOnParts) {
			setNumRow(nRow);
			pRow = permRow();
			pCol = permCol();
		}
		else {
			pRow = m_pPermutSparse[0].elementPntr();
			pCol = m_pPermutSparse[1].elementPntr();
		}
	}
	else {
		pRow = m_pPermutSparse[0].elementPntr();
		pCol = pPermCol;
	}

	setStabilizerLength(nRow - 1);
	setStabilizerLengthAut(ELEMENT_MAX);

	const auto len = nRow * sizeof(*pRow);
	memcpy(*pPermRows = pRow, m_pTrivialPermutCol, len);
	memcpy(pOrbits, m_pTrivialPermutCol, len);
	if (!pPermCol) {
		memcpy(pCol, m_pTrivialPermutCol, numCol() * sizeof(*pCol));

		if (!groupOnParts) {
			for (auto iPart = numParts; iPart--;) {
				auto pColPermStorage = permStorage(iPart);
				pColPermStorage->initPermutStorage();
				if (savePerm)
					pColPermStorage->savePermut(numRow(), permRow());
			}

			setGroupOrder(1);
			if (permColStorage() && (savePerm || permRowStorage())) {
				permColStorage()->initPermutStorage();
				if (savePerm)
					permColStorage()->savePermut(numCol(), pCol);
			}


			if (permRowStorage())
				permRowStorage()->initPermutStorage();
		}
	}

	return pCol;
}
CanonicityChecker(void)::addAutomorphism(const T numRow, const T *permRow, T *pOrbits, bool rowPermut, bool savePermut, bool calcGroupOrder)
{
	UpdateOrbits(permRow, numRow, pOrbits, rowPermut, calcGroupOrder);
	if (!rowPermut) {
		if (permRowStorage())
			permRowStorage()->savePermut(numRow, permRow);
	}
	else {
		if (savePermut)
			permStorage()->savePermut(numRow, permRow);
	}
}

#if USE_ASM <= 1   // We are not using Assembly OR we are using inline Assembly
CanonicityChecker(int)::checkColOrbit(T orbLen, T nColCurr, const S *pRow, const T *pRowPerm, T *permColumns) const
{
#if USE_ASM == 1
	_asm {	
			mov ecx, this
			// Reset counter
			mov esi, [ecx].m_pCounter			// counter()->element()
			mov eax, [esi + 4]	// ecx + 8		// numElement() - which is 2 for Insidence systems
			mov esi, [esi]		/// esi + 4
			bt eax, 0
			jz a1
			mov[esi + 4 * eax - 4], 0
       a1:  shr eax, 1
		    jmp a3
	   a2:  mov[esi + 8 * eax], 0
			mov[esi + 8 * eax + 4], 0
	   a3:  dec eax
			jge a2

			// Reset colNumbStorage()'s
			mov edx, [ecx].m_nColNumbStorage
			mov ebx, [ecx].m_rank				// i
			dec ebx								// since m_rank is always >= 2							
	   b1:  mov esi, [edx + ebx * 4]			// colNumbStorage()[i]
			mov [esi + 8], 0					// m_nNumb = 0
	        dec ebx
			jge b1

			// Count all entrances of the orbit:
			mov edi, [ecx].m_pPermutCol			// permCol()
			mov eax, [nColCurr]
			mov esi, [pRow]
			add esi, eax						// pRow += nColCurr
			shl eax, 2
			mov edi, [edi]	 // [edi + 4]
			add edi, eax						// permColumnCurr = permCol() + nColCurr
			mov ebx, [orbLen]
			dec ebx                             // since orbLen > 1
	   c1:  mov eax, [edi + 4*ebx]				// colNum = *(permColumnCurr + i)
			mov ecx, [pRowPerm]					// *(pRowPerm + colNum)
			add ecx, eax
			movsx eax, [ecx]					
			mov ecx, [edx + 4*eax]				//	colNumbStorage()[*(pRowPerm + colNum)]
			mov eax, [ecx].m_nNumb
			inc [ecx].m_nNumb  					// m_nNumb++
			mov ecx, [ecx] // ecx + 4                     // colNumbStorage()[*(pRowPerm + colNum)]
			shl eax, 2
			add ecx, eax
			mov eax, [edi + 4*ebx]
			mov [ecx], eax						// colNumbStorage()[*(pRowPerm + colNum)]->addElement(colNum)
			movsx eax, [esi + ebx]					// load 1 byte *(pRow + i)
			mov ecx, this
			mov ecx, [ecx].m_pCounter			// counter()->element()
			mov ecx, [ecx] // ecx + 4
			inc [ecx + eax*4]
	        dec ebx
			jge  c1

			// Check counters to see if current colOrbit is not the same as it was before
			mov ecx, this
			mov ebx, [ecx].m_rank				// i
			mov ecx, [ecx].m_pCounter
			mov ecx, [ecx]  // ecx + 4
			dec ebx
	   d1:  mov eax, [ecx + ebx*4]				// counter()->element(i)
		    mov esi, [edx + ebx*4]
			sub eax, [esi].m_nNumb
			je  d2
			mov eax, 1							// the orbit is lexicographycally larger 
			jg  done
	     	or eax, 0FFFFFFFFh                  // the orbit is lexicographycally smaller 
			jmp done
	   d2:  dec ebx
			jg  d1

			// Reordering of columns
			mov ecx, this
			mov ebx, [ecx].m_rank				// j
			dec ebx
	   e1:  mov ecx, [edx + ebx * 4]			// pStorage = colNumbStorage()[j];
			mov eax, [ecx].m_nNumb
			mov ecx, [ecx]						// pStorage->element
			jmp e3
	   e2:	mov esi, [ecx + eax * 4]            // pStorage->element(k)
			mov [edi], esi
			add edi, 4							// to next permCol element
	   e3:  dec eax
			jge e2
	        dec ebx
			jge e1
			xor eax, eax
     done:
	}
#else
	counter()->resetArray();

	for (auto i = rank(); i--;)
		colNumbStorage()[i]->resetArray();

	const auto *permColumnCurr = permColumns + nColCurr;
	pRow += nColCurr;
	for (auto i = orbLen; i--;) {
		const auto colNum = *(permColumnCurr + i);
		colNumbStorage()[*(pRowPerm + colNum)]->addElement(colNum);
		counter()->incCounter(*(pRow + i));
	}

	for (auto i = rank(); --i;) { // Don't need to check last element
		const auto diff = (int)counter()->element(i) - (int)colNumbStorage()[i]->numb();
		if (diff < 0)
			return -1;

		if (diff > 0)
			return 1;
	}

	// Reorder columns
	T i = nColCurr;
	for (auto j = rank(); j--;) {
		const auto pStorage = colNumbStorage()[j];
		for (auto k = pStorage->numb(); k--;)
			*(permColumns + i++) = pStorage->element(k);
	}

	return 0;
#endif
}
#endif

CanonicityChecker(void)::reconstructSolution(const ColOrbPntr pColOrbitStart, const ColOrbPntr pColOrbit,
	size_t colOrbLen, const ColOrbPntr pColOrbitIni, const T *pRowPerm, const T *pRowSolution, size_t solutionSize)
{
	// Skip all colOrbits which were equal to the tested solution
	int nOrb = 0;
	while (pColOrbitStart != pColOrbit) {
		nOrb++;
		pColOrbitStart = pColOrbitStart->next();
	}

	// Copying corresponding part of the solution:
	const auto rankM1 = rank() - 1;
	memcpy(improvedSolution(), pRowSolution, nOrb * rankM1 * sizeof(*improvedSolution()));
	auto *pImprovedSolution = improvedSolution() + nOrb * rankM1;
	memset(pImprovedSolution, 0, rankM1 * (solutionSize - nOrb) * sizeof(*pImprovedSolution));

	bool orbLenFlg = pColOrbit->length() > 1;
	size_t nColCurr = ((char *)pColOrbit - (char* )pColOrbitIni) / colOrbLen;
	while (true) {
		pImprovedSolution += rankM1;
		if (orbLenFlg) {
			// The elements in the counter make next (rank()-1) coordinates of our vector solution
			for (int i = rank(); --i;)
				*(pImprovedSolution - i) = static_cast<T>(colNumbStorage()[i]->numb());
		} else {
			// All next (rank()-1) coordinates except one of our vector solution are 0's. One, which is not 0 
			// is equal to 1 and is defined by corresponding element in the current row of the matrix
			const auto val = *(pRowPerm + *(permCol() + nColCurr));
			if (val > 0)
				*(pImprovedSolution - val) = 1;
		}

		pColOrbit = pColOrbit->next();
		if (!pColOrbit)
			return;

		const auto orbLen = pColOrbit->length();
		nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
		if ((orbLenFlg = (orbLen > 1))) {
			for (auto i = rank(); i--;)
				colNumbStorage()[i]->resetArray();

			// We need to calculate counters same way we did it in CCanonicityChecker::checkColOrbit()
			const auto *permColumnCurr = permCol() + nColCurr;
			for (auto i = orbLen; i--;) {
				const auto colNum = *(permColumnCurr + i);
				colNumbStorage()[*(pRowPerm + colNum)]->incNumb();
			}
		}
	}
}

CanonicityChecker(void)::outputAutomorphismInfo(FILE *file, const MatrixDataPntr pMatrix) const
{
	if (!(enumFlags() & t_outRowOrbits))
		return;         // We don't need the detailed group related information

	MUTEX_LOCK(out_mutex);
	const auto flag = enumFlags() & t_outRowPermute;
	outString(flag? "\nOrbits and generating permutations:\n" : "\nOrbits:\n", file);
	const auto pColOrbits = permColStorage() && (m_enumFlags & t_colOrbitsConstructed)? getColOrbits(0) : NULL;
	permStorage()->outputAutomorphismInfo(file, orbits(), permColStorage(), pColOrbits, pMatrix, flag);
	MUTEX_UNLOCK(out_mutex);
}

CanonicityChecker(bool)::groupIsTransitive() const
{ 
	auto nFixedRows = lenStabilizer();
	if (groupOrder() % (numRow() - nFixedRows))
		return false;

	const auto *pOrb = orbits();
	for (auto i = numRow(); i-- > nFixedRows;)
		if (*(pOrb + i) != nFixedRows)
			return false;

	return true;
}

CanonicityChecker(bool)::printMatrix(const designParam *pParam) const
{
#if PRINT_SOLUTIONS || PRINT_CURRENT_MATRIX
	if (printAll)
		return true;
#endif
	const auto outType = pParam->outType;
	return	outType & t_AllObject ||
			  outType & t_Transitive && groupIsTransitive() ||
			  outType & t_GroupOrderGT && groupOrder() > pParam->grpOrder ||
			  outType & t_GroupOrderLT && groupOrder() < pParam->grpOrder ||
			  outType & t_GroupOrderEQ && groupOrder() == pParam->grpOrder;
}

CanonicityChecker(T)::rowToChange(T nRow) const
{
	// Defines row of matrix which needs to be changed since it makes matrix non-canonical
	auto i = nRow;
    do {
		if (nRow < *(permRow() + i))
			nRow = *(permRow() + i);
    } while (i--);
    
	return nRow;
}

CanonicityChecker(T)::constructColIndex(const ColOrbPntr pColOrbit, const ColOrbPntr pColOrbitIni, size_t colOrbLen, T shift) const
{
	T idx = 0;
	while (pColOrbit) {
		// Define the number of columns to start with
		const auto numCol = shift + ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
		*(colIndex() + numCol) = idx++;
		pColOrbit = pColOrbit->next();
	}
	return idx;
}

size_t outString(const char *str, FILE *file)
{
    if (file)
        return fwrite(str, sizeof(*str), strlen(str), file) + 1;

	std::cout << str;
	return std::numeric_limits<std::size_t>::max();
}    

size_t outString(const char *str, const char *fileName, const char *mode)
{
	FOPEN(file, fileName, mode);
	const auto retVal = outString(str, file);
	FCLOSE(file);
	return retVal;
}
