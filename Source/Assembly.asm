	default rel
	global ?checkColOrbit@CCanonicityChecker@@AEAAHH_KPEBE1@Z

struc CPermut
	.m_pElem			resd	2	; pointer to T used un definition of CPermut as CSimpleArray<T>
	.m_nLen				resd	1	; integer member of CSimpleArray<T>
	.m_alignment		resd    1   ; we need it, because sizeof(CPermut) == 16
endstruc

struc CCanonicityChecker
	.CPermut			resd	4	; pointer to T used un definition of CPermut as CSimpleArray<T>
	.m_nStabLength		resd	1	; int 
    .m_rank				resd	1	; const int 
    .m_pPermutRow		resd	2	; CPermut *
    .m_pPermutCol		resd	2	; CPermut *
    .m_nColNumbStorage	resd	2	; CColNumbStorage **
	.m_pPermutStorage	resd	2	; CPermutStorage *
    .m_pCounter			resd	2	; CCounter<int> *
	.m_pColIndex		resd	2	; int *
	.m_pElem			resd	2	; pointer to T used un definition of CPermut as CSimpleArray<T>
	.m_nLen				resd	1	; integer member of CSimpleArray<T>
	.m_nGroupOrder:		resd	1	; uint 
	.m_nNumRow:			resd	1   ; int
	.m_nNumColOrb:		resd	1	; int
endstruc

struc CColNumbStorage
	.CPermut			resd    4
	.m_nNumb            resd    1
endstruc

;int CCanonicityChecker::checkColOrbit(int orbLen, size_t nColCurr, const MATRIX_ELEMENT_TYPE *pRow, const MATRIX_ELEMENT_TYPE *pRowPerm)
;{
;	counter()->resetArray();
;	for (int i = rank(); i--;)
;		colNumbStorage()[i]->resetArray();
;
;	const int *permColumnCurr = permCol() + nColCurr;
;	pRow += nColCurr;
;	for (int i = orbLen; i--;) {
;		const int colNum = *(permColumnCurr + i);
;		colNumbStorage()[*(pRowPerm + colNum)]->addElement(colNum);
;		counter()->incCounter(*(pRow + i));
;	}
;	for (int i = rank(); --i;) { // Don't need to check last element
;		const int diff = counter()->element(i) - colNumbStorage()[i]->numb();
;		if (diff < 0)
;			return -1;
;
;		if (diff > 0)
;			return 1;
;	}
;
;   // Reorder columns
;   size_t i = nColCurr;
;   for (int j = rank(); j--;) {
;       const CColNumbStorage *pStorage = colNumbStorage()[j];
;       for (int k = pStorage->numb(); k--;)
;			*(permCol() + i++) = pStorage->element(k);
;   }
;    
;	return 0;
;}

; 
; rsp+20h	- pRowPerm
; r9		- pRow
; r8		- nColCurr
; edx (rdx)	- orbLen
; rcx		- this

; To figured out the correct Assemler's name of the function,
; the technique from https://www.youtube.com/watch?v=UTNG-Ukvc3c was used
; Compile the code WITHOUT implementation of assembler's function,
; Linker will not find it, but will give the name it was expecting
; After that write the fubction with that name in assembler and do not forget to declare it as a "global" here

PROC_FRAME  ?checkColOrbit@CCanonicityChecker@@AEAAHH_KPEBE1@Z 
END_PROLOG
;			// Resert counter
			mov r11, [rcx + CCanonicityChecker.m_pCounter]	; counter()->element()
			movsxd rax, dword [r11 + CPermut.m_nLen]		; numElement()
			mov r11, [r11]
			bt eax, 0
			jz a1
			mov dword [r11 + 4 * rax - 4], 0
       a1:  shr rax, 1
		    jmp a3
	   a2:  mov long [r11 + 8 * rax], 0
	   a3:  dec rax
			jge a2

;			Reset colNumbStorage()'s
			mov r11, [rcx + CCanonicityChecker.m_nColNumbStorage]
			movsxd rax, dword [rcx + CCanonicityChecker.m_rank]	; i
			dec rax								; since m_rank is always >= 2							
	   b1:  mov r10, [r11 + 8*rax]				; colNumbStorage()[i]
			mov dword [r10 + CColNumbStorage.m_nNumb], 0	; m_nNumb = 0
	        dec rax
			jge b1

;			// Count all entrances of the orbit:
			mov r10, [rcx + CCanonicityChecker.m_pPermutCol] ; permCol()
			mov r10, [r10]
 			add r9, r8						    ; pRow += nColCurr
			shl r8, 2							;  4 * nColCurr
			add r10, r8						    ; permColumnCurr = permCol() + nColCurr
			mov r8, [rcx + CCanonicityChecker.m_pCounter]	; counter()->element()
			mov r8, [r8]
			dec rdx
	   c1:  movsxd rax, dword [r10 + 4*rdx]		; colNum = *(permColumnCurr + i)
			mov r11, [rsp + 20h + 8]			; since rsp was changed by 8 bytes when we call this function		
			add r11, rax						; *(pRowPerm + colNum)
			movsx eax, byte [r11]
			mov r11, [rcx + CCanonicityChecker.m_nColNumbStorage]							
			mov r11, [r11 + 8*rax]				;	colNumbStorage()[*(pRowPerm + colNum)]
			mov eax, [r11 + CColNumbStorage.m_nNumb]
			inc dword [r11 + CColNumbStorage.m_nNumb]	; m_nNumb++
			mov r11, [r11]                      ; colNumbStorage()[*(pRowPerm + colNum)]
			shl eax, 2
			add r11, rax
			mov eax, [r10 + 4*rdx]
			mov [r11], eax						; colNumbStorage()[*(pRowPerm + colNum)]->addElement(colNum)
			movsx eax, byte [r9 + rdx]			; load 1 byte *(pRow + i)
			inc dword [r8 + rax*4]
	   c2:  dec rdx
			jge  c1

;			Check counters to see if current colOrbit is not the same as it was before

			mov r11, [rcx + CCanonicityChecker.m_nColNumbStorage]
			mov r9, [rcx + CCanonicityChecker.m_pCounter]
			mov r9, [r9]
			movsxd rdx, dword [rcx + CCanonicityChecker.m_rank]	; i
			dec rdx
 	   d1:  mov eax, [r9 + rdx*4]				; counter()->element(i)
		    mov r8, [r11 + rdx*8]
			sub eax, [r8 + CColNumbStorage.m_nNumb]
 			je  d2
 			mov rax, 1							; the orbit is lexicographycally larger 
 			jg  done
 	     	or rax, 0FFFFFFFFh					; the orbit is lexicographycally smaller 
 			ret
 	   d2:  dec rdx
 			jg  d1

;			Reordering of columns
			movsxd rdx, dword [rcx + CCanonicityChecker.m_rank]	; j
            dec rdx
			mov r8, rdx
	   e1:  mov r9, [r11 + r8 * 8]						; pStorage = colNumbStorage()[j];
			mov edx, [r9 + CColNumbStorage.m_nNumb]
			mov r9, [r9]								; pStorage->element
			jmp e3
	   e2:	mov eax, [r9 + rdx * 4]						; pStorage->element(k)
			mov [r10], eax
			add r10, 4									; to next permCol element
	   e3:  dec rdx
			jge e2
	        dec r8
			jge e1
			xor rax, rax
     done:	ret
ENDPROC_FRAME
