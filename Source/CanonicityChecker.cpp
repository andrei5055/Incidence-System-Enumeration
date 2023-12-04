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

CanonicityChecker(T)::next_permutation(T *perm, const T *pOrbits, T idx, T lenStab) {
	// Function generates next permutation among those which stabilize first lenStab elements
	// We are using the algorithm from http://nayuki.eigenstate.org/res/next-lexicographical-permutation-algorithm/nextperm.java
	// taking into account that we don't need the the permutations which are equivalent with respect to already found orbits of the 
	// automorphis group acting on the matrix's rows.
	//
	// For instance, we found the automorphism
	//  (0, 1, 2, ..., i-1, pi, ...) and pi != i
	// it means that the elements pi and i are in the same orbits with respect to any stabilizer of the length < i.
	// It also means that on any place j 
	//  (0, 1, 2, ... , j-1,  j, ...) we don't need to try 
	//  (0, 1, 2, ... , j-1, pi, ...) when we already checked the permutation
	//  (0, 1, 2, ... , j-1,  i, ...)

	// Find non-increasing suffix
    const auto nRow = numRow();
	T temp, i, j;

	// Check if the algorithm, used immediately after 
	// some non-trivial automorphism was found
	const auto IDX_MAX = ELEMENT_MAX - 1;
	if (idx == IDX_MAX && perm[stabilizerLength()] == nRow - 1)
		idx = ELEMENT_MAX;

    if (idx == IDX_MAX) {
        // Firts call after some automorphism was found
        temp = perm[idx = i = stabilizerLength()];
        for (j = nRow; --j > temp;)
            perm[j] = j;
        
        for (auto k = j++; k-- > i;)
            perm[k+1] = k;
    } else {
        if (idx >= IDX_MAX) {
            j = i = nRow;
            while (--i > 0 && perm[i - 1] >= perm[i]);
        
            if (i == lenStab)
                return ELEMENT_MAX;

            // Find successor to pivot
            temp = perm[--i];
            while (perm[--j] <= temp);
        } else {
            temp = perm[j = i = idx];
            while (++j < nRow && perm[j] <= temp);
            if (j >= nRow) {
                revert(perm, nRow, i);
                return next_permutation(perm, pOrbits);
            }
        }
    }

    if (stabilizerLength() == i) {
        bool flag = false;
		auto k = j, tmp = perm[j];
        if (idx >= IDX_MAX) {
            while (k > i && *(pOrbits + perm[k]) != perm[k])
                k--;
            
            if (k != j) {
                if (!k)
                    return ELEMENT_MAX;
            
                flag = k == i;
                tmp = perm[k--];
                while (++k < j)
                    perm[k] = perm[k + 1];
            }
        } else {
            while (k < nRow && *(pOrbits + perm[k]) != perm[k])
                k++;
            
            if (k != j) {
                flag = k == nRow;
                if (flag) {
                    if (!i)
                        return ELEMENT_MAX;
                    
                    // Re-establish trivial permutation
                    k = idx - 1;
                    while (++k < j)
                        perm[k] = k;
                } else {
                    tmp = perm[k++];
                    while (--k > j)
                        perm[k] = perm[k - 1];
                }
            }
        }
        
        perm[j] = tmp;
        if (flag) {
            j = idx >= ELEMENT_MAX - 1? nRow - 1 : i;
            temp = perm[--i];
            setStabilizerLength(i);
        }
    }
    
    perm[i] = perm[j];
    perm[j] = temp;
    if (idx >= ELEMENT_MAX - 1) {
 		if (stabilizerLength() > i)
			setStabilizerLength(i);

		revert(perm, nRow, i);
	}

	return i;
}

CanonicityChecker(void)::UpdateOrbits(const T *permut, const T lenPerm, T *pOrb, bool rowPermut, bool calcGroupOrder)
{
	const T idx = udpdateStabLength(permut, lenPerm, pOrb, calcGroupOrder, rowPermut);
	permStorage()->UpdateOrbits(permut, lenPerm, pOrb, idx);
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
	S i = nColCurr;
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
	S nFixedRows = lenStabilizer();
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

CanonicityChecker(T)::nextPermutation(T *perm, const T *pOrbits, T idx, T lenStab) {
	// Function generates next permutation for the k-system 
	// Find non-increasing suffix
	const auto nRow = numRow();
	T temp, i, j;

	// Check if the algorithm, used immediately after 
	// some non-trivial automorphism was found
	const auto IDX_MAX = ELEMENT_MAX - 1;
	if (idx == IDX_MAX && perm[stabilizerLength()] == nRow - 1)
		idx = ELEMENT_MAX;

	if (idx == IDX_MAX) {
		// Firts call after some automorphism was found
		temp = perm[idx = i = stabilizerLength()];
		for (j = nRow; --j > temp;)
			perm[j] = j;

		for (auto k = j++; k-- > i;)
			perm[k + 1] = k;
	}
	else {
		if (idx >= IDX_MAX) {
			j = i = nRow;
			while (--i > 0 && perm[i - 1] >= perm[i]);

			if (i == lenStab)
				return ELEMENT_MAX;

			// Find successor to pivot
			temp = perm[--i];
			while (perm[--j] <= temp);
		}
		else {
			temp = perm[j = i = idx];
			while (++j < nRow && perm[j] <= temp);
			if (j >= nRow) {
				revert(perm, nRow, i);
				return nextPermutation(perm, pOrbits);
			}
		}
	}

	if (stabilizerLength() == i) {
		bool flag = false;
		auto k = j, tmp = perm[j];
		if (idx >= IDX_MAX) {
			while (k > i && *(pOrbits + perm[k]) != perm[k])
				k--;

			if (k != j) {
				if (!k)
					return ELEMENT_MAX;

				flag = k == i;
				tmp = perm[k--];
				while (++k < j)
					perm[k] = perm[k + 1];
			}
		}
		else {
			while (k < nRow && *(pOrbits + perm[k]) != perm[k])
				k++;

			if (k != j) {
				flag = k == nRow;
				if (flag) {
					if (!i)
						return ELEMENT_MAX;

					// Re-establish trivial permutation
					k = idx - 1;
					while (++k < j)
						perm[k] = k;
				}
				else {
					tmp = perm[k++];
					while (--k > j)
						perm[k] = perm[k - 1];
				}
			}
		}

		perm[j] = tmp;
		if (flag) {
			j = idx >= ELEMENT_MAX - 1 ? nRow - 1 : i;
			temp = perm[--i];
			setStabilizerLength(i);
		}
	}

	perm[i] = perm[j];
	perm[j] = temp;
	if (idx >= ELEMENT_MAX - 1) {
		if (stabilizerLength() > i)
			setStabilizerLength(i);

		revert(perm, nRow, i);
	}

	return i;
}


CanonicityChecker(bool)::rollBack(T *p_dayRes, T *p_dayIsUsed, int &j, int nDays) const
{
	while (j > 1) { // Do we need to go to previous day?
		p_dayIsUsed[p_dayRes[--j]] = 0;
		if (++p_dayRes[j] < nDays) {
			j--;
			return true;
		}

		p_dayRes[j] = 0;
	}

	return false;
}

CanonicityChecker(bool)::CheckCanonicity(const T *result, int nDays) {
	// return true  - continue, 
	//        false - stop (calculate new matrix)
	const auto lenStab = stabiliserLengthExt();
/*
	const auto* pEnum = pCanonParam->pEnum;
	auto* pPartNumb = pCanonParam->pPartNumb;
	auto* pRowOut = pCanonParam->pRowOut;
	const auto* pGroupOnParts = pCanonParam->pGroupOnParts;
	// Construct trivial permutations for rows and columns
	const auto rowPermut = outInfo & t_saveRowPermutations;
	const auto* pMatr = pCanonParam->pMatrix;
	const auto* pMatrPerm = pMatr;
	const auto numParts = pCanonParam->numParts;
	bool check_trivial_row_perm = pCanonParam->pPermCol != NULL;
	auto savePermut = rowPermut && (enumFlags() & t_outRowPermute);

	const auto colOrbLen = pEnum->colOrbitLen();
	const auto* pPartInfo = pMatr->partsInfo();
	const auto nonCombinedDesign = numParts == 1;

	// Reset permutation counters, if used
	PREPARE_PERM_OUT(permColStorage());
	PREPARE_PERM_OUT(permRowStorage());

	bool calcGroupOrder = true;
	bool retVal = true;
	bool usingGroupOnBlocks = false;

	size_t numGroups = 0;
	T startingRowNumb = pCanonParam->startingRowNumb;
	size_t idxPerm[16] = {};	// to keep the indices of currently used permutation
	size_t* pIndxPerms = NULL;	// of i-th symmetrical group acting on the parts

	T idxPartSrc[16] = {};		// to keep the initial (source) indices of the parts
	T* pPartSrc = numParts <= countof(idxPartSrc) ? idxPartSrc : new T[numParts];
	for (auto i = numParts; i--;)
		pPartSrc[i] = i;

	const S* pCurrSolution;
	size_t solutionSize;
	if (pRowSolution) {
		// Because the position of current solution (pRowSolution->solutionIndex()) 
		// can be changed, take pointer here
		pCurrSolution = pRowSolution->currSolution();
		solutionSize = pRowSolution->solutionLength();
#if USE_STRONG_CANONICITY
		solutionStorage()->clear();
		pRowSolution->setLenOrbitOfSolution(0);
#endif
	}
*/
	static int cntr = 0;
	size_t startIndex = 0;
	T* permColumn = NULL;
	auto pOrbits = orbits();
	const auto* res = result;
	const auto lenMemory = m_numElem + 2 * nDays;
	T buff[100];
	T *p_players = lenMemory <= countof(buff) ? buff : new T[lenMemory];
	T *p_dayRes = p_players + m_numElem;
	T *p_dayIsUsed = p_dayRes + nDays;

	const auto lenGroup = rank();
	const auto numGroup = m_numElem / lenGroup;
	for (int iDay = 0; iDay < nDays; iDay++, res += lenGroup) {
		if (res[0] || !copyTuple(res, p_players))
			return false;

		T inc = 0;
		for (auto j = numGroup; --j;) {
			if (!copyTuple(res += lenGroup, p_players, inc+=lenGroup) ||
				p_players[*res] < p_players[*(res - lenGroup)]) // Comparing first elements of the groups
				return false;
		}

		if (lenGroup >= 3) {
			// Ordering last two elements of each 
			auto *res = p_players;
			for (auto j = numGroup; j--; res += lenGroup) {
				if (res[1] > res[2]) {  // TODO: Write code for lenGroup > 3
					const auto tmp = res[1];
					res[1] = res[2];
					res[2] = tmp;
				}
			}
		}

		// Check canonicity of the codes for the other days
		memset(p_dayRes, 0, 2 * nDays * sizeof(*p_dayRes));
		p_dayIsUsed[p_dayRes[0] = iDay] = 1;
        int j = 0;
		T k;
		while (true) {
			while (++j < nDays) {
				cntr++;
				// Looking for the first unused day
				k = p_dayRes[j];
				while (k < nDays && p_dayIsUsed[k])
					k++;

				if (k == nDays) // Unable to find unused day.
					break;

				p_dayIsUsed[p_dayRes[j] = k] = 1;
				const auto *resDayPerm = result + k * m_numElem;
				const auto *resDay = result + j * m_numElem;
				int diff = 0;
				T t = -1;
				while (++t < m_numElem && !(diff = (int)p_players[resDayPerm[t]] - resDay[t]))
					;
				if (t < m_numElem) {
					if (diff < 0)
						return false;
					else
						break;
				}
			}

			if (rollBack(p_dayRes, p_dayIsUsed, j, nDays))
				continue;
//			else

/*
			if (j > 1) { // Do we need to go to previous day?
				p_dayIsUsed[p_dayRes[--j]] = 0;
				p_dayRes[j--]++;
				continue;
			}
*/
			if (true || k == nDays) {// Unable to find unused day.
				break;   // there is no previous day for which a different choice can be made 
			}
			// automorphism found
		}

		p_dayIsUsed[iDay] = 0;
		continue;   // temporary
		T* permPlayers = NULL;
		permColumn = init(m_numElem, 0, false, pOrbits, &permPlayers, false, permColumn);

		T nElem = ELEMENT_MAX;
/*
		if (check_trivial_row_perm) {
			nRow = startingRowNumb;
			goto try_permut;
		}
		*/
		while (true) {

//		next_permut:
			nElem = nextPermutation(permPlayers, pOrbits, nElem, lenStab);
			if (nElem == ELEMENT_MAX || nElem < lenStabilizer())
				break;

			for (T iDay = 0; iDay < nDays;  iDay++) {
				const auto* pDayRes = result + iDay * m_numElem;
				for (; nElem < m_numElem; nElem++) {
					// const auto* pRow = pMatr->GetRow(nElem);
				}
			}
		}
	}

	if (p_players != buff)
		delete[] p_players;

	return true;
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
