
#include "InsSysEnumerator.h"

#if USE_THREADS
	#include <iostream>
	#include <system_error>
	#if USE_BOOST
		#include <boost/exception_ptr.hpp>
		// How to install the C++ Boost Libraries on Windows
		// http://andres.jaimes.net/718/how-to-install-the-c-boost-libraries-on-windows/

		#if USE_POOL
			// Using advices from: stackoverflow.com/questions/19500404/how-to-create-a-thread-pool-using-boost-in-c
			#include <boost/asio/io_service.hpp>
			#include <boost/thread/thread_pool.hpp>
		#endif
	#else
		#include <thread>
	#endif

	using namespace THREAD_NAME_SPACE;
#endif

#include "InSysSolver.h"
#include "EnumInfo.h"
#include "GroupsInfo.h"
#ifdef _MSC_VER
#ifndef _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC
#endif
#include <crtdbg.h>
#endif
#if defined(_MSC_VER) && defined(_DEBUG)
#define new new(_NORMAL_BLOCK, THIS_FILE, __LINE__)
#undef THIS_FILE
	static char THIS_FILE[] = __FILE__;
#endif
#if PRINT_SOLUTIONS || PRINT_CURRENT_MATRIX
int ccc = 0;
#endif

#if USE_MUTEX
std::mutex out_mutex;
#endif

CEnumerator::CEnumerator(const CMatrix *pMatrix, bool matrOwner) : m_pMatrix(pMatrix), CColOrbitManager(pMatrix->maxElement() + 1, pMatrix->rowNumb(), pMatrix->colNumb())
{ 
	const auto rowNumb = pMatrix->rowNumb();
	m_pRow = new CRowSolution *[rowNumb];
	setRowEquation(NULL);
	setCanonChecker(new CCanonicityChecker(rowNumb, colNumb(), rank()));
    setOutFile(NULL);
	setMatrOwner(matrOwner);
}

CEnumerator::~CEnumerator()
{
	delete[] rowStuff(rowMaster());
	delete[] rowStuffPntr();
	delete rowEquation();
    delete canonChecker();
	if (isMatrOwner())
		delete matrix();
}

CRowSolution *CEnumerator::FindRowSolution(PERMUT_ELEMENT_TYPE lastRightPartIndex)
{
	const size_t nVar = MakeSystem();
	if (nVar == -1)
        return NULL;
    
	CRowSolution *pNextRowSolution = FindSolution(nVar, lastRightPartIndex);
    if (!pNextRowSolution)
        return NULL;

    OUTPUT_SOLUTION(pNextRowSolution, outFile(), false);
	return sortSolutions(pNextRowSolution, lastRightPartIndex) ? pNextRowSolution : NULL;
}

#if USE_THREADS
#if WRITE_MULTITHREAD_LOG
void thread_message(int threadIdx, const char *pComment, t_threadCode code, void *pntr = 0)
{
	FILE *file = fopen("Report.txt", "a");
	fprintf(file, "thrIdx = %2d: %10s  (%d) %p\n", threadIdx, pComment, code, pntr);
	fclose(file);
}
#endif

void threadEnumerate(CThreadEnumerator *threadEnum, const CEnumerator *pMaster)
{
	threadEnum->EnumerateBIBD(pMaster);
}

int CEnumerator::threadWaitingLoop(int thrIdx, t_threadCode code, CThreadEnumerator *threadEnum, int nThread) const
{
	// Try to find the index of not-running thread
	size_t loopingTime = 0;
	while (true) {
		if (loopingTime > REPORT_INTERVAL) {
			// We run this loop enough to send report message
			enumInfo()->reportProgress(threadEnum, nThread);
			loopingTime = 0;
		}

		const int startIdx = thrIdx;
#if WRITE_MULTITHREAD_LOG 
		FILE *file = fopen("Report.txt", "a");
		fprintf(file, "startIdx = %2d  Spipping code = %d\n", startIdx, code);
#endif
		while (threadEnum[thrIdx].code() == code) {
			if (++thrIdx == nThread)
				thrIdx = 0;

			if (thrIdx == startIdx)
				break;
		}

		CThreadEnumerator *pEnum = threadEnum + thrIdx;
#if WRITE_MULTITHREAD_LOG 
		fprintf(file, "==>thrIdx = %2d  CODE = %d\n", thrIdx, pEnum->code());
		fclose(file);
#endif
		switch (pEnum->code()) {
			case t_threadFinished:
				enumInfo()->updateGroupInfo(pEnum->enumInfo());
				thread_message(thrIdx, "finished", code);
				enumInfo()->reportProgress(pEnum);
				if (code == t_threadNotUsed) {
					pEnum->setCode(t_threadNotUsed);
					continue;
				}

			case t_threadLaunchFail:
				pEnum->reInit();

			case t_threadNotUsed:
				thread_message(thrIdx, "notUsed", code);
				return thrIdx;

			case t_threadRunning:
				if (++thrIdx == nThread)
					thrIdx = 0;

			case t_threadUndefined:
				loopingTime += SLIP_TIME;
				this_thread::sleep_for(chrono::microseconds(SLIP_TIME));
				break;
		}
	}
}
#endif

ulonglong CEnumerator::Enumerate(bool writeFile, CEnumInfo *pEnumInfo, const CEnumerator *pMaster, t_threadCode *pTreadCode)
{
	char buff[256];
	if (writeFile) {
		if (makeFileName(buff, countof(buff))) {
			FOPEN(file, buff, "w");
			setOutFile(file);
			makeJobTitle(buff, countof(buff), "\n");
			outString(buff, outFile());
		}
        
		if (pEnumInfo && !pMaster && makeFileName(buff, countof(buff), "tmp.txt"))
			pEnumInfo->setReportFileName(buff);
	}

	CEnumInfo enumInfoInst;
	if (!pEnumInfo)
		pEnumInfo = &enumInfoInst;

#if SOLUTION_STATISTICS
	int nSol, nCntr, nFirst;
	uint nMax;
	nSol = nCntr = nFirst = nMax = 0;
#endif

#if USE_THREADS
	CThreadEnumerator *pThreadEnum = NULL;
	int thrIdx = 0;
	#if USE_POOL
		asio::io_service *pIoService = NULL;
		thread_group *pThreadpool = NULL;
	#endif
#endif

	// Allocate memory for the orbits of two consecutive rows
	size_t nRow;
	CRowSolution *pRowSolution;
	InitRowSolutions(pMaster);
	const bool threadFlag = pMaster != NULL;
	if (threadFlag) {
		setCurrentRowNumb(nRow = pMaster->currentRowNumb());
		pRowSolution = pMaster->rowStuff(nRow);
		setOutFile(pMaster->outFile());
		setX0_3(pMaster->getX0_3());
		const auto firstUnforced = pMaster->firstUnforcedRow();
		if (firstUnforced > 0) {
			setFirstUnforcedRow(firstUnforced);
			memcpy(forcibleLambdaPntr() + firstUnforced, pMaster->forcibleLambdaPntr() + firstUnforced, (rowNumb() - firstUnforced) * sizeof(*forcibleLambdaPntr()));
		}
	} else {
		nRow = 0;
		pRowSolution = setFirstRowSolutions();
		setEnumInfo(pEnumInfo);
		pEnumInfo->startClock();

#if USE_THREADS 
		pThreadEnum = new CThreadEnumerator[USE_THREADS];
	#if USE_POOL
		// Create an asio::io_service and a thread_group (through pool in essence)
		pIoService = new asio::io_service();
		pThreadpool = new thread_group();

		// This will start the ioService processing loop.All tasks
		// assigned with ioService.post() will start executing.
		asio::io_service::work work(*pIoService);

		// This will add USE_THREADS threads to the thread pool.
		for (auto i = 0; i < USE_THREADS; i++)
			(pThreadEnum + i)->setThread(pThreadpool->create_thread(bind(&asio::io_service::run, pIoService)));

//		pThreadpool->join_all();
	#endif
#endif
	}
    
	if (!threadFlag)
		prepareToTestExtraFeatures();

	// For multithreaded version we need to test only one top level solution
	const size_t nRowEnd = nRow ? nRow + 1 : 0;
    initiateColOrbits(rowNumb(), isIS_enumerator(), pMaster);
	PERMUT_ELEMENT_TYPE level;
	while (pRowSolution) {

		const bool useCanonGroup = USE_CANON_GROUP && nRow > 0;

#if USE_THREADS
		if (!nRowEnd && nRow == MT_LEVEL) {
#if WRITE_MULTITHREAD_LOG
			for (int i = USE_THREADS; i--;) threadEnum[i].setThreadID(i);
#endif
			// We are in master enumerator
			while (pRowSolution) {
				(pThreadEnum+thrIdx)->setupThreadForBIBD(this, nRow);
				do {
					try {
#if USE_POOL
						pIoService->post(bind(threadEnumerate, pThreadEnum + thrIdx, this));
						pIoService->run_one();
						pIoService->stop();
//						(pThreadEnum + thrIdx)->getThread()->detach();
//						pThreadpool->join_all();
#else
						thread t1(threadEnumerate, pThreadEnum + thrIdx, this);
						t1.detach();
#endif
						thread_message(thrIdx, "detached", (pThreadEnum + thrIdx)->code());
						break;
					}
					catch (const std::system_error& e) {
						std::cout << "System_error with code " << e.code()
							<< " meaning " << e.what() << " thrIdx = " << thrIdx << '\n';
					}
#if USE_BOOST
					catch (exception const& ex) {
						exception_ptr except = current_exception();
						char const** pInfo = exception_detail::get_info< throw_function >::get(ex);
						std::cout << "Boost exception" << std::endl;
						std::cout << "==>" << pInfo[1] << std::endl;
					}
#endif
					// So far I saw only one type of exception: 
					// e.code() = "generic:11 (EAGAIN)
					// e.what() = "meaning resource unavailable try again :"
					// Let's wait. Perhaps, it will help.
					this_thread::sleep_for(chrono::seconds(1));
				} while (true);

				// Try to find the index of not-running thread
				thrIdx = threadWaitingLoop(thrIdx, t_threadRunning, pThreadEnum, USE_THREADS);
				pRowSolution = pRowSolution->NextSolution(useCanonGroup);
			}

#if WAIT_THREADS
			// All canonocal solutions are distributed amongst the threads
			// Waiting for all threads finish their jobs
			threadWaitingLoop(thrIdx, t_threadNotUsed, pThreadEnum, USE_THREADS);
			thread_message(999, "DONE", t_threadUndefined);
			pEnumInfo->reportProgress(t_reportNow);
#else
			pEnumInfo->reportProgress(pThreadEnum, USE_THREADS);
#endif
		} else {
#endif
#if PRINT_SOLUTIONS
			ccc++;
#endif		
			OUTPUT_SOLUTION(pRowSolution, outFile(), true);
			const VECTOR_ELEMENT_TYPE *pCurrSolution = pRowSolution->currSolution();
			CColOrbit *pColOrb = MakeRow(pCurrSolution);
			if (nRow == 2)
				setX0_3(*pCurrSolution);

			OUTPUT_MATRIX(matrix(), outFile(), nRow + 1);
			if (++nRow == rowNumb()) {
				pEnumInfo->incrConstrTotal();
				bool flag = true;
				if (canonChecker()->TestCanonicity(nRow, this, t_saveRowToChange + t_saveRowPermutations, &level)) {
//					Construct Aut(D)
//					int ddd = canonChecker()->constructGroup();
					if (TestFeatures(pEnumInfo)) {
						if (noReplicatedBlocks() && pEnumInfo->constructedAllNoReplBlockMatrix()) {
							pEnumInfo->setNoReplBlockFlag(false);
							level = getInSys()->GetK();
							flag = false;
						}
						else {
							pEnumInfo->updateConstrCounters(this);
							if (PRINT_MATRICES)
								matrix()->printOut(outFile(), nRow, pEnumInfo->constrCanonical(), canonChecker());

#if PRINT_SOLUTIONS
							ccc = 0;
#endif
							if (!rowMaster())  // We are not in the slave thread
								pEnumInfo->reportProgress(t_matrConstructed);
						}
					}
				}
				else
					flag = false;

				if (!flag) {
					while (--nRow > level) {
						setCurrentRowNumb(nRow);
						rowStuff(nRow)->resetSolution();
						resetColOrbitCurr();
						resetUnforcedColOrb();
					}
					pRowSolution = rowStuff(nRow);
					setCurrentRowNumb(nRow);
				}
				else {
					nRow--;
					pRowSolution = NULL;
				}
			}
			else {
				setCurrentRowNumb(nRow);
				setColOrbitCurr(pColOrb);
				setCurrUnforcedOrbPtr(nRow);
				if (!USE_CANON_GROUP || canonChecker()->TestCanonicity(nRow, this, 0, &level, pRowSolution)) {
					// We need to get lastRightPartIndex here and use later because 
					// for multi-thread configuration it could be changed by master
					const PERMUT_ELEMENT_TYPE lastRightPartIndex = pRowSolution->solutionIndex();
					if (pMaster) {
						copyInfoFromMaster(pMaster);
#if WAIT_THREADS
						pMaster = NULL;
						*pTreadCode = t_threadRunning;
#endif
					}

					if (!useCanonGroup)
						canonChecker()->setGroupOrder(1);

					pRowSolution = FindRowSolution(lastRightPartIndex);
#if USE_THREADS && !WAIT_THREADS
					if (pMaster && pRowSolution) {
						pMaster = NULL;
						*pTreadCode = t_threadRunning;
					}
#endif
					checkUnusedSolutions(pRowSolution);

#if SOLUTION_STATISTICS
					if (pRowSolution) {
						nCntr++;
						nSol += pRowSolution->numSolutions();
						nFirst += lastRightPartIndex;
						if (nMax < pRowSolution->numSolutions())
							nMax = pRowSolution->numSolutions();
					}
#endif

					OUTPUT_CANON_GROUP(useCanonGroup, canonChecker(), outFile());
				} else {
					if (pMaster)
						break;

					pRowSolution = NULL;
				}
			}
#if USE_THREADS
		}
#endif
		while (!pRowSolution || !(pRowSolution = pRowSolution->NextSolution(useCanonGroup))) {
            rowStuff(nRow)->resetSolution();
            resetColOrbitCurr();
            resetUnforcedColOrb();
			if (nRow-- > nRowEnd) {
				setCurrentRowNumb(nRow);
				pRowSolution = rowStuff(nRow);
			} else
				break;
		}
	} 

#if USE_THREADS && !WAIT_THREADS
	if (!threadFlag)
		threadWaitingLoop(thrIdx, t_threadNotUsed, pThreadEnum, USE_THREADS);
#endif

    closeColOrbits();

	if (!threadFlag || !USE_THREADS) { // We are not in the slave thread
		pEnumInfo->outEnumInfo(outFile());

#if SOLUTION_STATISTICS
		SPRINTF(buff, "nCntr = %5d:  %10.1f - %10.1f  nMax = %5d\n", nCntr, ((float)nSol) / nCntr, ((float)nFirst) / nCntr, nMax);
		outString(buff, outFile());
#endif
	} else {
		if (pTreadCode)
			*pTreadCode = pMaster ? t_threadLaunchFail : t_threadFinished;
	}

#if USE_THREADS
	delete[] pThreadEnum;
#endif
	return pEnumInfo->constrCanonical();
} 

CColOrbit *CEnumerator::MakeRow(const VECTOR_ELEMENT_TYPE *pRowSolution) const
{
    const auto nRow = currentRowNumb();
	const bool nextColOrbNeeded = nRow + 1 < rowNumb();
	auto *pRow = matrix()->GetRow(nRow);
	memset(pRow, 0, sizeof(*pRow) * colNumb());

	const CColOrbit *pColOrbit = colOrbit(nRow);
    CColOrbit *pNextRowColOrbit = colOrbit(nRow+1);
    
	const int maxElement = rank() - 1;
    const size_t colOrbLen = colOrbitLen();
	const CColOrbit *pColOrbitIni = colOrbitIni(nRow);
    CColOrbit *pNextRowColOrbitNew = NULL;
    CColOrbit *pColOrbitLast = NULL;
	while (pColOrbit) {
        CColOrbit *pNewColOrbit = NULL;
		// Define the number of column to start with
		const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
        auto lenRemaining = pColOrbit->lenght();
		auto *pRowCurr = pRow + nColCurr;
		for (int i = rank(); i--;) {
			const MATRIX_ELEMENT_TYPE lenFragm = i? pRowSolution[maxElement - i] : lenRemaining;
			if (!lenFragm)
				continue;
            
			if (nextColOrbNeeded) {
				if (!pNewColOrbit) {
					pNewColOrbit = (CColOrbit *)((char *)pNextRowColOrbit + nColCurr * colOrbLen);
					if (pColOrbitLast)
						pColOrbitLast->setNext(pNewColOrbit);
					else
						pNextRowColOrbitNew = pNewColOrbit;
				}

				pNewColOrbit = (pColOrbitLast = pNewColOrbit)->InitOrbit(lenFragm, colOrbLen, pColOrbit, i);
			}

            if (!i)
                break;
            
            // Construct corresponding part of matrix's nRow's row
#if MATRIX_ELEMENT_IS_BYTE
			memset(pRowCurr, i, lenFragm);
            pRowCurr += lenFragm;
#else
			for (int j = 0; j < lenFragm; j++)
				*pRowCurr[j] = i;
#endif
            if (!(lenRemaining -= lenFragm))
                break;    // remainin parts of the current solution are 0's
		}
        
		pRowSolution += maxElement;
		pColOrbit = pColOrbit->next();
	}

	if (getUnforceColOrbPntr()) {
        // Set unforced elements:
        for (auto row = firstUnforcedRow(); row <= nRow; row++) {
			const CColOrbit *pColOrbitIni = colOrbitIni(row);
			CColOrbit **ppUnforced = unforcedOrbits(row);
            for (int i = 1; i < rank(); i++) {
                pColOrbit = *(ppUnforced + i);
                while (pColOrbit) {
                    const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
                    const auto lenFragm = pColOrbit->lenght();
                    if (lenFragm > 1) {
#if MATRIX_ELEMENT_IS_BYTE
                        memset(pRow + nColCurr, i, lenFragm);
#else
                        for (int j = 0; j < lenFragm; j++)
                            pRow[nColCurr + j] = i;
#endif
                    } else
                        pRow[nColCurr] = i;
                    
                    pColOrbit = pColOrbit->next();
                }
            }
        }
    }
    
    if (pColOrbitLast)    // Not the last row of the matrix
        pColOrbitLast->setNext(NULL);
    
    return pNextRowColOrbitNew;
}

void CEnumerator::InitRowSolutions(const CEnumerator *pMaster)
{
	const auto nRow = pMaster? pMaster->currentRowNumb() + 1 : 0;
	CRowSolution *pSolutions = new CRowSolution[rowNumb() - nRow];
	for (auto i = rowNumb(); i-- > nRow;)
		m_pRow[i] = pSolutions + i - nRow;

	if (nRow) 
		memcpy(rowStuffPntr(), pMaster->rowStuffPntr(), nRow * sizeof(*rowStuffPntr()));
}

#if USE_STRONG_CANONICITY_A
void CEnumerator::checkUnusedSolutions(CRowSolution *pRowSolution)
{
	if (!pRowSolution)
		return;

	// We just about to test first solution for row nRow. Before we do that
	// we could test solutions which SHOULD NOT be used for this row (because with their usages
	// the matrix could not be constructed to the last row). But these solutions 
	//  a) could be eliminated since the do not satisfy the conditions of strong canonicity; 
	//  b) they could be moved to some new orbits of solutions
	const PERMUT_ELEMENT_TYPE solIdx = pRowSolution->solutionIndex();
	if (solIdx == -1)
		return;

	const int nRow = currentRowNumb() + 1;
	pRowSolution->setSolutionIndex(-1);
	CRowSolution *pRowSol = pRowSolution;
	size_t level = solIdx + 1;
	while (pRowSol->NextSolution(true) && pRowSol->solutionIndex() <= solIdx) {
		const VECTOR_ELEMENT_TYPE *pCurrSol = pRowSol->currSolution();
		CColOrbit *pColOrb = MakeRow(pCurrSol);

		setCurrentRowNumb(nRow);
		setCurrUnforcedOrbPtr(nRow);
		setColOrbitCurr(pColOrb);
		canonChecker()->TestCanonicity(nRow, this, 0, &level, pRowSol);
		rowStuff(nRow)->resetSolution();
		resetColOrbitCurr();
		resetUnforcedColOrb();
		setCurrentRowNumb(nRow - 1);
	}

	pRowSolution->setSolutionIndex(level - 1);
}
#endif
