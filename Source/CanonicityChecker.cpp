//
//  CanonicityChecker.cpp
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 3/28/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#include "Enumerator.h"

template class CCanonicityChecker<TDATA_TYPES>;

CanonicityChecker(void)::InitCanonicityChecker(S nRow, S nCol, int rank, S *pMem)
{
	m_rank = rank;
	m_pPermutRow = (CPermut *)(pMem += nRow);
	m_pPermutRow->Init(nRow, pMem = (S *)((char *)pMem + sizeof(CPermut)));
	m_pPermutCol = (CPermut *)(pMem += nRow);
	m_pPermutCol->Init(nCol, pMem = (S *)((char *)pMem + sizeof(CPermut)));
	setColIndex(pMem += nCol);
	m_pCounter = (CCounter<int> *)(pMem += (nCol << 1));
	m_pCounter->Init(rank, (int *)(pMem = (S *)((char *)pMem + sizeof(CCounter<int>))));
	setPermStorage((PermutStoragePntr)(pMem = (S *)((char *)pMem + rank * sizeof(int))));
	m_nColNumbStorage = (CColNumbStorage **)(pMem = (S *)((char *)pMem + sizeof(Class2(CPermutStorage))));
	pMem = (S *)((char *)pMem + rank * sizeof(CColNumbStorage *));
	for (int i = rank; i--;) {
		m_nColNumbStorage[i] = (CColNumbStorage *)(pMem);
		pMem = (S *)((char *)pMem + sizeof(CColNumbStorage));
		m_nColNumbStorage[i]->Init(nCol, pMem);
		pMem += nCol;
	}

	m_pImprovedSol = NULL;
	setSolutionStorage(NULL);
}

CanonicityChecker(void)::init(S nRow, bool savePerm) {
	setNumRow(nRow);
	setStabilizerLength(nRow - 1);
	setStabilizerLengthAut(ELEMENT_MAX);

	auto* pRow = permRow();
	for (auto i = nRow; i--;)
		*(pRow + i) = *(orbits() + i) = i;

	auto* pCol = permCol();
	for (auto nCol = numCol(); nCol-- > nRow;)
		*(pCol + nCol) = nCol;

	memcpy(pCol, pRow, nRow * sizeof(*permCol()));
	for (auto iPart = numParts(); iPart--;) {
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

CanonicityChecker(void)::revert(S i)
{
    // Reverse suffix (if needed)
    auto *array = permRow();
    auto j = numRow();
    while (++i < --j)
        array[i] ^= (array[j] ^= (array[i] ^= array[j]));
}

CanonicityChecker(S)::next_permutation(S idx, S lenStab) {
	// Funstion generates next permutation amoungth those which stabilize first lenStab elements
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
    auto *array = permRow();
    const auto nRow = numRow();
	S temp, i, j;

	// Check if the algorithm, used immediately after 
	// some automorphism was found
	const auto IDX_MAX = ELEMENT_MAX - 1;
	if (idx == IDX_MAX && array[stabilizerLength()] == nRow - 1)
		idx = ELEMENT_MAX;

    if (idx == IDX_MAX) {
        // Firts call after some automorphism was found
        temp = array[idx = (int)(i = stabilizerLength())];
        for (j = nRow; --j > temp;)
            array[j] = j;
        
        for (auto k = j++; k-- > i;)
            array[k+1] = k;
    } else {
        if (idx >= IDX_MAX) {
            j = i = nRow;
            while (--i > 0 && array[i - 1] >= array[i]);
        
            if (i == lenStab)
                return ELEMENT_MAX;

            // Find successor to pivot
            temp = array[--i];
            while (array[--j] <= temp);
        } else {
            temp = array[j = i = idx];
            while (++j < nRow && array[j] <= temp);
            if (j >= nRow) {
                revert(i);
                return next_permutation();
            }
        }
    }

    if (stabilizerLength() == i) {
        bool flag = false;
		auto k = j, tmp = array[j];
        if (idx >= IDX_MAX) {
            while (k > i && *(orbits() + array[k]) != array[k])
                k--;
            
            if (k != j) {
                if (!k)
                    return ELEMENT_MAX;
            
                flag = k == i;
                tmp = array[k--];
                while (++k < j)
                    array[k] = array[k + 1];
            }
        } else {
            while (k < nRow && *(orbits() + array[k]) != array[k])
                k++;
            
            if (k != j) {
                flag = k == nRow;
                if (flag) {
                    if (!i)
                        return ELEMENT_MAX;
                    
                    // Re-establish trivial permutation
                    k = idx - 1;
                    while (++k < j)
                        array[k] = k;
                } else {
                    tmp = array[k++];
                    while (--k > j)
                        array[k] = array[k - 1];
                }
            }
        }
        
        array[j] = tmp;
        if (flag) {
            j = idx >= ELEMENT_MAX - 1? nRow - 1 : i;
            temp = array[--i];
            setStabilizerLength(i);
        }
    }
    
    array[i] = array[j];
    array[j] = temp;
    if (idx >= ELEMENT_MAX - 1) {
 		if (stabilizerLength() > i)
			setStabilizerLength(i);

		revert(i);
	}

	return i;
}

CanonicityChecker(void)::UpdateOrbits(const S *permut, S lenPerm, S *pOrb, bool rowPermut, bool calcGroupOrder)
{
	S idx = 0;
	while (idx == permut[idx])
		idx++;

	if (calcGroupOrder) {
		if (rowPermut && stabilizerLengthAut() > idx)
			updateGroupOrder();

		setStabilizerLength(idx);
		setStabilizerLengthAut(idx);
	}

	permStorage()->UpdateOrbits(permut, lenPerm, pOrb, idx);
}

CanonicityChecker(void)::addAutomorphism(bool rowPermut, bool savePermut)
{
	if (!rowPermut) {
		if (permRowStorage())
			permRowStorage()->savePermut(numRow(), permRow());
	}
	else {
		UpdateOrbits(permRow(), numRow(), orbits(), rowPermut, true);
		if (savePermut)
			permStorage()->savePermut(numRow(), permRow());
	}
}

CanonicityChecker(void)::updateGroupOrder()
{
	S len = 1;
	const auto *pOrb = orbits();
	const auto i = stabilizerLengthAut();
	auto idx = i;
	while (++idx < numRow()) {
		if (*(pOrb + idx) == i)
			len++;
	}

	setGroupOrder(len * groupOrder());
}

#if USE_ASM <= 1   // We are not using Assembly OR we are using inline Assembly
CanonicityChecker(int)::checkColOrbit(S orbLen, S nColCurr, const T *pRow, const T *pRowPerm) const
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

	auto *permColumns = permCol();
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
	size_t colOrbLen, const ColOrbPntr pColOrbitIni, const T *pRowPerm, const S *pRowSolution, size_t solutionSize)
{
	// Skip all colOrbits which were equal to the tested solution
	int nOrb = 0;
	while (pColOrbitStart != pColOrbit) {
		nOrb++;
		pColOrbitStart = pColOrbitStart->next();
	}

	// Copying corresponding part of the solution:
	const int rankM1 = rank() - 1;
	memcpy(improvedSolution(), pRowSolution, nOrb * rankM1 * sizeof(*improvedSolution()));
	auto *pImprovedSolution = improvedSolution() + nOrb * rankM1;
	memset(pImprovedSolution, 0, rankM1 * (solutionSize - nOrb) * sizeof(*pImprovedSolution));

	bool orbLenFlg = pColOrbit->length() > 1;
	size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
	while (true) {
		if (orbLenFlg) {
			// The elements in the counter make next (rank()-1) coordinates of our vector solution
			for (int i = rank(); --i;)
				*(pImprovedSolution + rankM1 - i) = static_cast<S>(colNumbStorage()[i]->numb());
		} else {
			// All next (rank()-1) coordinates except one of our vector solution are 0's. One, which is not 0 
			// is equal to 1 and is defined by corresponding element in the current row of the matrix
			const auto val = *(pRowPerm + *(permCol() + nColCurr));
			if (val > 0)
				*(pImprovedSolution + rankM1 - val) = 1;
		}

		pImprovedSolution += rankM1;
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
	if (groupOrder() % numRow())
		return false;

	const auto *pOrb = orbits();
	for (auto i = numRow(); i--;)
		if (*(pOrb + i))
			return false;

	return true;
}

CanonicityChecker(bool)::printMatrix(const designParam *pParam) const
{
	const auto outType = pParam->outType;
	return	outType & t_AllObject ||
			  outType & t_Transitive && groupIsTransitive() ||
			  outType & t_GroupOrderGT && groupOrder() > pParam->grpOrder ||
			  outType & t_GroupOrderLT && groupOrder() < pParam->grpOrder ||
			  outType & t_GroupOrderEQ && groupOrder() == pParam->grpOrder;
}

CanonicityChecker(S)::rowToChange(S nRow) const
{
	// Defines row of matrix which needs to be changed since it makes matrix non-canonical
	auto i = nRow;
    do {
		if (nRow < *(permRow() + i))
			nRow = *(permRow() + i);
    } while (i--);
    
	return nRow;
}

CanonicityChecker(S)::constructColIndex(const ColOrbPntr pColOrbit, const ColOrbPntr pColOrbitIni, size_t colOrbLen, S shift) const
{
	S idx = 0;
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