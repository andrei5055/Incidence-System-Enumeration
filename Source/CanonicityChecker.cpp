//
//  CanonicityChecker.cpp
//  BIBD_Mac
//
//  Created by Andrei Ivanov on 3/28/14.
//  Copyright (c) 2014 Andrei Ivanov. All rights reserved.
//

#include "Enumerator.h"

CCanonicityChecker::CCanonicityChecker(size_t nRow, size_t nCol, int rank) : CPermut(nRow), m_rank(rank)
{
    m_pPermutRow = new CPermut(nRow);
    m_pPermutCol = new CPermut(nCol);
	#pragma warning(suppress: 6297)
	setColIndex(new size_t [nCol<<1]);
    m_pCounter = new CCounter<int>(rank);
	setPermStorage(new CPermutStorage());
    m_nColNumbStorage = new CColNumbStorage *[rank];
    for (int i = rank; i--;)
        m_nColNumbStorage[i] = new CColNumbStorage(nCol);   

	m_pImprovedSol = new VECTOR_ELEMENT_TYPE[(rank - 1) * nCol];
	setSolutionStorage(new CSolutionStorage());
}

CCanonicityChecker::~CCanonicityChecker()
{
	delete m_pPermutRow;
	delete m_pPermutCol;
    delete counter();
	delete permStorage();
    for (int i = rank(); i--;)
        delete m_nColNumbStorage[i];
    
    delete [] colNumbStorage();
	delete[] colIndex();
	delete [] improvedSolution();
	delete solutionStorage();
}

void CCanonicityChecker::init(size_t nRow, bool savePerm)
{
	setNumRow(nRow);
    setStabilizerLength(nRow-1);
	setStabilizerLengthAut(-1);
    auto *pRow = permRow();
    for (size_t i = nRow; i--;)
        *(pRow + i) = *(orbits() + i) = i;
    
    auto *pCol = permCol();
    auto nCol = numCol();
    while (--nCol >= nRow)
        *(pCol + nCol) = nCol;
    
    memcpy(pCol, pRow, nRow * sizeof(*permCol()));
	initPermStorage();
	setGroupOrder(1);
	if (savePerm) 
		savePerm(numRow(), permRow());
}

void CCanonicityChecker::revert(size_t i)
{
    // Reverse suffix (if needed)
    auto *array = permRow();
    auto j = numRow();
    while (++i < --j)
        array[i] ^= (array[j] ^= (array[i] ^= array[j]));
}

size_t CCanonicityChecker::next_permutation(size_t idx) {
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
    size_t temp, i, j;

	// Check if the algorithm, used immediately after 
	// some automorphism was found
	if (idx == SIZE_MAX - 1 && array[stabilizerLength()] == nRow - 1)
		idx = SIZE_MAX;		// no, we will use the standart one

    if (idx == SIZE_MAX - 1) {
        // Firts call after some automorphism was found
        temp = array[idx = (int)(i = stabilizerLength())];
        for (j = nRow; --j > temp;)
            array[j] = j;
        
        for (auto k = j++; k-- > i;)
            array[k+1] = k;
    } else {
        if (idx >= SIZE_MAX - 1) {
            j = i = nRow;
            while (--i > 0 && array[i - 1] >= array[i]);
        
            if (i-- == 0)
                return -1;

            // Find successor to pivot
            temp = array[i];
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
        size_t k, tmp = array[k = j];
        if (idx >= SIZE_MAX - 1) {
            while (k > i && *(orbits() + array[k]) != array[k])
                k--;
            
            if (k != j) {
                if (!k)
                    return -1;
            
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
                        return -1;
                    
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
            j = idx >= SIZE_MAX - 1? nRow - 1 : i;
            temp = array[--i];
            setStabilizerLength(i);
        }
    }
    
    array[i] = array[j];
    array[j] = temp;
    if (idx >= SIZE_MAX - 1) {
 		if (stabilizerLength() > i)
			setStabilizerLength(i);

		revert(i);
	}

	return i;
}

void CCanonicityChecker::addAutomorphism(bool rowPermut)
{
    const auto *array = permRow();
    size_t idx = 0;
    while (idx == array[idx])
        idx++;

	if (rowPermut && stabilizerLengthAut() > idx)
		updateGroupOrder();

    setStabilizerLength(idx);
	setStabilizerLengthAut(idx);
    auto *pOrb = orbits();
    const auto nRow = numRow();
    do  {
        auto i = *(pOrb + idx);
        auto j = *(pOrb + array[idx]);
        if (j == i)
            continue;
        
        if (j < i) {
            i ^= j;
            i ^= j ^= i;
        }
        
        for (size_t k = 0; k < nRow; k++) {
            if (*(pOrb+k) == j)
                *(pOrb+k) = i;
        }
    } while (++idx < nRow);
    

	if (!rowPermut) {
		size_t *permCol;
		const auto lenPerm = getLenPermutCol(&permCol);
		if (groupOrder() == 1) {		// no permutations were saved yet
			savePerm(lenPerm, NULL);	// save trivial permutation
			setGroupOrder(2);			// 2 is a fake number, we use it just not be here again
		}

		savePerm(lenPerm, permCol);
	} else
		savePerm(numRow(), permRow());
}

void CCanonicityChecker::updateGroupOrder()
{
	const auto i = stabilizerLengthAut();
	auto idx = i;
	int len = 1;
	auto *pOrb = orbits();
	while (++idx < numRow()) {
		if (*(pOrb + idx) == i)
			len++;
	}

	setGroupOrder(len * groupOrder());
}

#if USE_ASM <= 1   // We are not using Assembly OR we are using inline Assembly
int CCanonicityChecker::checkColOrbit(size_t orbLen, size_t nColCurr, const MATRIX_ELEMENT_TYPE *pRow, const MATRIX_ELEMENT_TYPE *pRowPerm) const
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

	for (int i = rank(); i--;)
		colNumbStorage()[i]->resetArray();

	auto *permColumns = permCol();
	const auto *permColumnCurr = permColumns + nColCurr;
	pRow += nColCurr;
	for (auto i = orbLen; i--;) {
		const auto colNum = *(permColumnCurr + i);
		colNumbStorage()[*(pRowPerm + colNum)]->addElement(colNum);
		counter()->incCounter(*(pRow + i));
	}

	for (int i = rank(); --i;) { // Don't need to check last element
		const int diff = counter()->element(i) - colNumbStorage()[i]->numb();
		if (diff < 0)
			return -1;

		if (diff > 0)
			return 1;
	}

	// Reorder columns
	size_t i = nColCurr;
	for (int j = rank(); j--;) {
		const CColNumbStorage *pStorage = colNumbStorage()[j];
		for (int k = pStorage->numb(); k--;)
			*(permColumns + i++) = pStorage->element(k);
	}

	return 0;
#endif
}
#endif

void CCanonicityChecker::reconstructSolution(const CColOrbit *pColOrbitStart, const CColOrbit *pColOrbit, size_t colOrbLen, const CColOrbit *pColOrbitIni, const MATRIX_ELEMENT_TYPE *pRowPerm, const VECTOR_ELEMENT_TYPE *pRowSolution, size_t solutionSize)
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
	VECTOR_ELEMENT_TYPE *pImprovedSolution = improvedSolution() + nOrb * rankM1;
	memset(pImprovedSolution, 0, rankM1 * (solutionSize - nOrb) * sizeof(*pImprovedSolution));

	bool orbLenFlg = pColOrbit->lenght() > 1;
	size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
	while (true) {
		if (orbLenFlg) {
			// The elements in the counter make next (rank()-1) coordinates of our vector solution
			for (int i = rank(); --i;)
				*(pImprovedSolution + rankM1 - i) = colNumbStorage()[i]->numb();
		} else {
			// All next (rank()-1) coordinates except one of our vector solution are 0's. One, which is not 0 
			// is equal to 1 and is defined by corresponding element in the current row of the matrix
			const MATRIX_ELEMENT_TYPE val = *(pRowPerm + *(permCol() + nColCurr));
			if (val > 0)
				*(pImprovedSolution + rankM1 - val) = 1;
		}

		pImprovedSolution += rankM1;
		pColOrbit = pColOrbit->next();
		if (!pColOrbit)
			return;

		const auto orbLen = pColOrbit->lenght();
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

void CCanonicityChecker::outputAutomorphismInfo(FILE *file) const
{
	MUTEX_LOCK(out_mutex);
	outString("\nOrbits and generating permutations:\n", file);
	outputPerm(file, orbits(), numRow());
	outputPermutations(file, numRow());
	MUTEX_UNLOCK(out_mutex);
}

bool CCanonicityChecker::groupIsTransitive() const
{ 
	if (groupOrder() % numRow())
		return false;

	const size_t *pOrb = orbits();
	for (auto i = numRow(); i--;)
		if (*(pOrb + i))
			return false;

	return true;
}

bool CCanonicityChecker::TestCanonicity(size_t nRowMax, const CEnumerator *pEnum, int outInfo, size_t *pRowOut, CRowSolution *pRowSolution)
{
	// Construct trivial permutations for rows and columns
	const bool rowPermut = outInfo & t_saveRowPermutations;
	init(nRowMax--, rowPermut);

	const CMatrix *pMatr = pEnum->matrix();
	const auto colOrbLen = pEnum->colOrbitLen();
	const auto colNumb = pEnum->colNumb();
	CColOrbit **colOrbit = pEnum->colOrbits();
	CColOrbit **colOrbitIni = pEnum->colOrbitsIni();
	const VECTOR_ELEMENT_TYPE *pCurrSolution;
	size_t solutionSize;
	if (pRowSolution) {
		// Because the position of current solution (pRowSolution->solutionIndex()) 
		// could be changed, let's take pointer here
		pCurrSolution = pRowSolution->currSolution();
		solutionSize = pRowSolution->solutionSize();
#if USE_STRONG_CANONICITY
		solutionStorage()->clear();
		pRowSolution->setLenOrbitOfSolution(0);
#endif
	}

	const size_t *pColIndex = NULL;
	size_t *pVarPerm = NULL;
	bool retVal = true;
	size_t startIndex = 0;
	size_t nRow = SIZE_MAX;
	while (true) {
	next_permut:
		nRow = next_permutation(nRow);
		if (nRow == SIZE_MAX)
			break;

		OUTPUT_PERMUTATION(this, pEnum->outFile(), orbits(), numRow());

		// Loop for all remaining matrix's rows
		for (; nRow <= nRowMax; nRow++) {
			const MATRIX_ELEMENT_TYPE *pRow = pMatr->GetRow(nRow);
			const MATRIX_ELEMENT_TYPE *pRowPerm = pMatr->GetRow(*(permRow() + nRow));
			const CColOrbit *pColOrbitIni = colOrbitIni[nRow];
			const CColOrbit *pColOrbit = colOrbit[nRow];
 			while (pColOrbit) {

				// Define the number of column to start with
				const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
				const auto orbLen = pColOrbit->lenght();
				int diff;
				if (orbLen > 1)
					diff = checkColOrbit(orbLen, nColCurr, pRow, pRowPerm);
				else
					diff = (int)*(pRow + nColCurr) - *(pRowPerm + *(permCol() + nColCurr));

				if (diff < 0) {
#if RECURSIVE_CANON
					if (orbLen > 1)
						return false;

					// Construct matrix, for just found permRow() and permCol()
					CMatrix *pMatr = matrix();
					CMatrix matr(pMatr, permRow(), permCol());
					CCanonicityChecker canon(rowNumb(), colNumb(), rank());
					setMatrix(&matr);
					std::cout << "########";
					matr.printOut();
					setCanonChecker(&canon);

					// Check, if this matrix is canonical
					int level;
					bool retVal = TestCanonicity(level);
					setMatrix(pMatr);
#else
					if (pRowOut) {
						if (outInfo & t_saveRowToChange)
							rowToChange(pRowOut, nRow);
						else {
							if (pRowSolution && !pRowSolution->isLastSolution()) {
								// NOTE: No need to remove last solution, we will leave this level anyway
								if (nRow < nRowMax) {
									// We can remove all solutions which are in the same group with the solution just tested.
									// We don't need them even as the right parts of our equations, because the usage of everyone 
									// will make the matrix non canonical
									pRowSolution->removeNoncanonicalSolutions(startIndex);
								} else {
									reconstructSolution(colOrbit[nRow], pColOrbit, colOrbLen, pColOrbitIni, pRowPerm, pCurrSolution, solutionSize);
#if USE_STRONG_CANONICITY
									if (solutionStorage()) {
										for (auto *pSolution : *solutionStorage()) {
											if (!memcmp(pSolution, improvedSolution(), solutionSize*sizeof(*pSolution)))
												goto next_permut; // We already saw this solution
										}
									}

									startIndex = pRowSolution->moveNoncanonicalSolutions(improvedSolution(), startIndex, solutionStorage(), pRowOut);
									if (startIndex != SIZE_MAX) {
										retVal = false;
										// we will try to find better alternative
										goto next_permut;
									}
#else
									pRowSolution->moveNoncanonicalSolutions(improvedSolution(), 0, solutionStorage());
#endif
								}
							}
						}
					}

					return false;
#endif
				}

				if (diff > 0)
					goto next_permut;

				pColOrbit = pColOrbit->next();
			}
		}

		// Automorphism found:
#if (PRINT_CURRENT_MATRIX && PRINT_PERMUTATION)
		outString("-----\n", pEnum->outFile());
#endif
		if (!rowPermut) {
			// We are here to define the canonicity of partially constructed matrix AND
			// we just found the nontrivial automorphism.
			const CColOrbit *pColOrbitIni = colOrbitIni[nRow];
			const CColOrbit *pColOrbit = colOrbit[nRow];
			if (!pColIndex) // Index for columns was not yet constructed
				pVarPerm = (size_t *)((pColIndex = constructColIndex(pColOrbit, pColOrbitIni, colOrbLen)) + colNumb);

			size_t varIdx = 0;
			while (pColOrbit) {
				const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
				*(pVarPerm + varIdx++) = *(pColIndex + *(permCol() + nColCurr));
				pColOrbit = pColOrbit->next();
			}
		}

		addAutomorphism(rowPermut);
		nRow = SIZE_MAX - 1;
	}

	if (rowPermut)
		updateGroupOrder();

	return retVal;
}

void CCanonicityChecker::rowToChange(size_t *pRowOut, size_t nRow) const
{
	// Defines row of matrix which needs to be changed since it makes matrix non-canonical
    size_t i = nRow;
    do {
		if (nRow < *(permRow() + i))
			nRow = *(permRow() + i);
    } while (i--);
    
	*pRowOut = nRow;
}

size_t *CCanonicityChecker::constructColIndex(const CColOrbit *pColOrbit, const CColOrbit *pColOrbitIni, size_t colOrbLen)
{
	int idx = 0;
	while (pColOrbit) {
		// Define the number of column to start with
		const size_t numCol = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
		*(colIndex() + numCol) = idx++;
		pColOrbit = pColOrbit->next();
	}

	setNumColOrb(idx);
	return colIndex();
}

size_t CCanonicityChecker::getLenPermutCol(size_t **permCol)
{
	*permCol = colIndex() + numCol();
	return numColOrb();
}

void outString(const char *str, FILE *file)
{
    if (file)
        fputs(str, file);
    else
		std::cout << str;
}    

void outString(const char *str, const char *fileName, const char *mode)
{
	FOPEN(file, fileName, mode);
	outString(str, file);
	FCLOSE(file);
}