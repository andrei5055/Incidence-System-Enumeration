
#include "InsSysEnumerator.h"
#include "CanonicityChecker.h"
#include "ThreadEnumerator.h"

#if USE_THREADS_ENUM
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
#include <sys/stat.h>
#if defined(_WIN32) || defined(_WIN64)
	#include <direct.h>
	#define MKDIR(x) _mkdir(x)
#elif defined(__linux__)
	#define MKDIR(x)  mkdir(dirName), 0777)
#else
	"Please define corresponding method for making directory"
#endif

#include <mutex>  

#if PRINT_SOLUTIONS || PRINT_CURRENT_MATRIX
int ccc = 0;
#endif

#if USE_MUTEX
std::mutex out_mutex;
#endif

bool fileExists(const char *path, bool file = true)
{
	struct stat info;
	return (stat(path, &info) == 0) && (file || info.st_mode & S_IFDIR);
}

template class CEnumerator<MATRIX_ELEMENT_TYPE>;

template<class T>
CRowSolution<T> *CEnumerator<T>::FindRowSolution(PERMUT_ELEMENT_TYPE lastRightPartIndex)
{
	prepareToFindRowSolution();
	const size_t nVar = MakeSystem();
	if (nVar == (size_t )-1)
        return NULL;
    
	setPrintResultNumVar(nVar);
	auto pNextRowSolution = FindSolution(nVar, lastRightPartIndex);
    if (!pNextRowSolution)
        return NULL;

    OUTPUT_SOLUTION(pNextRowSolution, outFile(), false);
	return sortSolutions(pNextRowSolution, lastRightPartIndex) ? pNextRowSolution : NULL;
}

#if USE_THREADS_ENUM
#if WRITE_MULTITHREAD_LOG
void thread_message(int threadIdx, const char *pComment, t_threadCode code, void *pntr)
{
	FILE *file = fopen("Report.txt", "a");
	fprintf(file, "thrIdx = %2d: %10s  (%d) %p\n", threadIdx, pComment, code, pntr);
	fclose(file);
}
#endif

template<typename T>
void threadEnumerate(CThreadEnumerator<T> *threadEnum, designRaram *param, const CEnumerator<T> *pMaster)
{
	threadEnum->EnumerateBIBD(param, pMaster);
}

template<class T>
int CEnumerator<T>::threadWaitingLoop(int thrIdx, t_threadCode code, CThreadEnumerator<T> *threadEnum, size_t nThread) const
{
	// Try to find the index of not-running thread
	size_t loopingTime = 0;
	while (true) {
		if (loopingTime > REPORT_INTERVAL) {
			// We run this loop enough to send report message
			this->enumInfo()->reportProgress(threadEnum, nThread);
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

		CThreadEnumerator<T> *pEnum = threadEnum + thrIdx;
#if WRITE_MULTITHREAD_LOG 
		fprintf(file, "==>thrIdx = %2d  CODE = %d\n", thrIdx, pEnum->code());
		fclose(file);
#endif
		switch (pEnum->code()) {
			case t_threadFinished:
				LAUNCH_CANONICITY_TESTING(pEnum, this);
				this->enumInfo()->updateGroupInfo(pEnum->enumInfo());
				thread_message(thrIdx, "finished", code);
				this->enumInfo()->reportProgress(pEnum);
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

template<class T>
ulonglong CEnumerator<T>::Enumerate(designRaram *pParam, bool writeFile, CEnumInfo<T> *pEnumInfo, const CEnumerator<T> *pMaster, t_threadCode *pTreadCode)
{
#if !CONSTR_ON_GPU
	std::mutex mtx;
	char buff[256], jobTitle[256];
	const size_t lenBuffer = countof(buff);
	const int mt_level = pParam->mt_level;
	const size_t threadNumb = pParam->threadNumb;

	size_t lenName = 0;
	if (writeFile) {
		// We will be here only for the master
		pParam->firstMatr = true;
		if (!makeJobTitle(jobTitle, countof(jobTitle), "\n"))
			return (size_t)-1;

		// Construct the file name of the file with the enumeration results
		if (!makeFileName(buff, lenBuffer, ""))
			return (size_t)-1;

		lenName = strlen(buff);

		// Create file name for the output of final enumeration results
		strcpy_s(buff + lenName, lenBuffer - lenName, FILE_NAME(""));

		// The results are known, if the file with the enumeration results exists
		const bool knownResults = fileExists(buff);
		if (knownResults)
			strcpy_s(buff + lenName, lenBuffer - lenName, FILE_NAME(CURRENT_RESULTS));
		else
			lenName = 0;

		// Create a new file for output of the enumeration results
		FOPEN(file, buff, "w");
		this->setOutFile(file);
		outString(jobTitle, this->outFile());
	}

	// Create file name for the output of intermediate results 
	if ((writeFile || pEnumInfo && !pMaster) && makeFileName(buff, lenBuffer, FILE_NAME(INTERMEDIATE_RESULTS)))
		pEnumInfo->setReportFileName(buff);
#endif

#if SOLUTION_STATISTICS
	int nSol, nCntr, nFirst;
	uint nMax;
	nSol = nCntr = nFirst = nMax = 0;
#endif

#if USE_THREADS_ENUM
	CThreadEnumerator<T> *pThreadEnum = NULL;
	int thrIdx = 0;
	#if USE_POOL
		asio::io_service *pIoService = NULL;
		thread_group *pThreadpool = NULL;
	#endif
#endif
	auto pMatrix = static_cast<const CMatrix<T> *>(this->matrix());

	// Allocate memory for the orbits of two consecutive rows
	T nRow;
	CRowSolution<T> *pRowSolution;
	InitRowSolutions(pMaster);
	const bool threadFlag = pMaster != NULL;
	if (threadFlag) {
		this->setCurrentRowNumb(nRow = pMaster->currentRowNumb());
		pRowSolution = pMaster->rowStuff(nRow);
		this->setOutFile(pMaster->outFile());
		setX0_3(pMaster->getX0_3());
		const auto firstUnforced = pMaster->firstUnforcedRow();
		if (firstUnforced > 0) {
			setFirstUnforcedRow(firstUnforced);
			memcpy(forcibleLambdaPntr() + firstUnforced, pMaster->forcibleLambdaPntr() + firstUnforced, (rowNumb() - firstUnforced) * sizeof(*forcibleLambdaPntr()));
		}
	} else {
		nRow = 0;
		pRowSolution = setFirstRowSolutions();
		this->setEnumInfo(pEnumInfo);
		pEnumInfo->startClock();

#if USE_THREADS_ENUM 
		pThreadEnum = new CThreadEnumerator<T>[pParam->threadNumb];
	#if USE_POOL
		// Create an asio::io_service and a thread_group (through pool in essence)
		pIoService = new asio::io_service();
		pThreadpool = new thread_group();

		// This will start the ioService processing loop.All tasks
		// assigned with ioService.post() will start executing.
		asio::io_service::work work(*pIoService);

		// This will add threadNumb threads to the thread pool.
		for (auto i = 0; i < threadNumb; i++)
			(pThreadEnum + i)->setThread(pThreadpool->create_thread(bind(&asio::io_service::run, pIoService)));

//		pThreadpool->join_all();
	#endif
#endif
	}
    
	if (!threadFlag)
		prepareToTestExtraFeatures();

	// For multi-threaded version we need to test only one top level solution
	const size_t nRowEnd = nRow ? nRow + 1 : 0;
    this->initiateColOrbits(rowNumb(), this->IS_enumerator(), pMaster);
	T level;
	while (pRowSolution) {

		const bool useCanonGroup = USE_CANON_GROUP && nRow > 0;

#if USE_THREADS_ENUM
		if (!nRowEnd && nRow == mt_level) {
#if WRITE_MULTITHREAD_LOG
			for (int i = threadNumb; i--;) threadEnum[i].setThreadID(i);
#endif
			// We are in master enumerator
			while (pRowSolution) {
				(pThreadEnum+thrIdx)->setupThreadForBIBD(this, nRow, thrIdx);
				do {
					try {
#if USE_POOL
						pIoService->post(bind(threadEnumerate, pThreadEnum + thrIdx, pParam, this));
						pIoService->run_one();
						pIoService->stop();
//						(pThreadEnum + thrIdx)->getThread()->detach();
//						pThreadpool->join_all();
#else
						thread t1(threadEnumerate<T>, pThreadEnum + thrIdx, pParam, this);
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
				thrIdx = threadWaitingLoop(thrIdx, t_threadRunning, pThreadEnum, threadNumb);
				pRowSolution = pRowSolution->NextSolution(useCanonGroup);
			}

#if WAIT_THREADS
			// All canonocal solutions are distributed amongst the threads
			// Waiting for all threads finish their jobs
			threadWaitingLoop(thrIdx, t_threadNotUsed, pThreadEnum, threadNumb);
			thread_message(999, "DONE", t_threadUndefined);
			pEnumInfo->reportProgress(t_reportNow);
#else
			pEnumInfo->reportProgress(pThreadEnum, threadNumb);
#endif
		} else {
#endif
#if PRINT_SOLUTIONS
			ccc++;
#endif		
			REPORT_PROGRESS(pEnumInfo, t_reportByTime);
			OUTPUT_SOLUTION(pRowSolution, outFile(), true);
			const auto *pCurrSolution = pRowSolution->currSolution();
			auto *pColOrb = MakeRow(pCurrSolution);
			if (nRow == 2)
				setX0_3(*pCurrSolution);

			OUTPUT_MATRIX(pMatrix, outFile(), nRow + 1);
			if (++nRow == rowNumb()) {
				pEnumInfo->incrConstrTotal();
				bool flag = true;
				if (!TestCanonicityOnGPU()) {
					EXIT(-1);
					if (this->TestCanonicity(nRow, this, t_saveRowToChange + t_saveRowPermutations, &level)) {
						//					Construct Aut(D)
						//					int ddd = canonChecker()->constructGroup();
						int matrFlags;
						if (TestFeatures(pEnumInfo, this->matrix(), &matrFlags)) {
							if (noReplicatedBlocks() && pEnumInfo->constructedAllNoReplBlockMatrix()) {
								pEnumInfo->setNoReplBlockFlag(false);
								level = getInSys()->GetK();
								flag = false;
							}
							else {
								pEnumInfo->updateConstrCounters(matrFlags, this->groupOrder(), this->groupIsTransitive());
#if !CONSTR_ON_GPU
								if (this->printMatrix(pParam)) {
									mtx.lock();
									if (pParam->firstMatr) {
										pParam->firstMatr = false;
										outString(BEG_OUT_BLOCK "Constructed Matrices: " END_OUT_BLOCK, this->outFile());
									}

									pMatrix->printOut(this->outFile(), nRow, pEnumInfo->constrCanonical(), this);
									mtx.unlock();
								}
#endif
#if PRINT_SOLUTIONS
								ccc = 0;
#endif
								if (!this->rowMaster())  // We are not in the slave thread
									REPORT_PROGRESS(pEnumInfo, t_matrConstructed);
							}
						}
					}
					else
						flag = false;
				}

				if (!flag) {
					while (--nRow > level) {
						this->setCurrentRowNumb(nRow);
						rowStuff(nRow)->resetSolution();
						this->resetColOrbitCurr();
						this->resetUnforcedColOrb();
					}
					pRowSolution = rowStuff(nRow);
					this->setCurrentRowNumb(nRow);
				}
				else {
					nRow--;
					pRowSolution = NULL;
				}
			}
			else {
				this->setCurrentRowNumb(nRow);
				this->setColOrbitCurr(pColOrb);
				this->setCurrUnforcedOrbPtr(nRow);
				if (!USE_CANON_GROUP || this->TestCanonicity(nRow, this, 0, &level, pRowSolution)) {
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
						this->setGroupOrder(1);

					setPrintResultRowNumber(nRow);
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
#if USE_THREADS_ENUM
		}
#endif
		while (!pRowSolution || !(pRowSolution = pRowSolution->NextSolution(useCanonGroup))) {
            rowStuff(nRow)->resetSolution();
            this->resetColOrbitCurr();
            this->resetUnforcedColOrb();
			if (nRow-- > nRowEnd) {
				this->setCurrentRowNumb(nRow);
				pRowSolution = rowStuff(nRow);
			} else
				break;
		}
	} 

#if USE_THREADS_ENUM && !WAIT_THREADS
	if (!threadFlag)
		threadWaitingLoop(thrIdx, t_threadNotUsed, pThreadEnum, pParam->threadNumb);
#endif

    this->closeColOrbits();

	const ulonglong retVal = pEnumInfo->constrCanonical();
	if (!threadFlag || !USE_THREADS_ENUM) {

#if USE_THREADS_ENUM
		delete[] pThreadEnum;
#else
//#ifdef USE_CUDA
		// This method is called after thread is ended, When they are used
		LAUNCH_CANONICITY_TESTING(enumInfo(), this);
//#endif
#endif

#if !CONSTR_ON_GPU
		// We are not in the slave thread
		if (!pParam->firstMatr)
			outString(END_OUT_BLOCK "Constructed Matrices " BEG_OUT_BLOCK "\n", this->outFile());

		pEnumInfo->outEnumInfo(this->outFilePntr(), lenName == 0);

		t_resType resType;
		if (lenName) {
			// Compare current results with previously obtained
			bool betterResults = true;
			const char *currentFile = FILE_NAME(CURRENT_RESULTS);
			strcpy_s(buff + lenName, countof(buff) - lenName, currentFile);
			if (compareResults(buff, lenName, &betterResults)) {
				// Create the file name with the current results 
				strcpy_s(jobTitle, buff);
				strcpy_s(jobTitle + lenName, countof(jobTitle) - lenName, currentFile);

				if (betterResults) {
					remove(buff);			// Remove file with previous results
					rename(jobTitle, buff);	// Rename file
					resType = t_resBetter;
				}
				else {
					if (pParam->firstMatr)
						remove(jobTitle);	// Deleting new file only when it does not contain matrices

					resType = t_resWorse;
				}
			}
			else
				resType = t_resInconsistent; // results are not the same as before
		}
		else
			resType = t_resNew;

		pEnumInfo->setResType(resType);

		if (pParam->outType & t_Summary)
			pEnumInfo->outEnumInformation(this->outFilePntr());

		this->setEnumInfo(NULL);
#endif

#if SOLUTION_STATISTICS
		SPRINTF(buff, "nCntr = %5d:  %10.1f - %10.1f  nMax = %5d\n", nCntr, ((float)nSol) / nCntr, ((float)nFirst) / nCntr, nMax);
		outString(buff, outFile());
#endif
	} else {
		if (pTreadCode)
			*pTreadCode = pMaster ? t_threadLaunchFail : t_threadFinished;
	}

	return retVal;
} 

template<class T>
CColOrbit<T> *CEnumerator<T>::MakeRow(const VECTOR_ELEMENT_TYPE *pRowSolution) const
{
    const auto nRow = this->currentRowNumb();
	const bool nextColOrbNeeded = nRow + 1 < rowNumb();
	auto *pRow = this->matrix()->GetRow(nRow);
	memset(pRow, 0, sizeof(*pRow) * this->colNumb());

	const auto *pColOrbit = this->colOrbit(nRow);
    auto *pNextRowColOrbit = this->colOrbit(nRow+1);
    
	const int maxElement = this->rank() - 1;
    const auto colOrbLen = this->colOrbitLen();
	const auto *pColOrbitIni = this->colOrbitIni(nRow);
    CColOrbit<T> *pNextRowColOrbitNew = NULL;
    CColOrbit<T> *pColOrbitLast = NULL;
	while (pColOrbit) {
        CColOrbit<T> *pNewColOrbit = NULL;
		// Define the number of column to start with
		const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
        auto lenRemaining = pColOrbit->length();
		auto *pRowCurr = pRow + nColCurr;
		for (int i = this->rank(); i--;) {
			const auto lenFragm = i? pRowSolution[maxElement - i] : lenRemaining;
			if (!lenFragm)
				continue;
            
			if (nextColOrbNeeded) {
				if (!pNewColOrbit) {
					pNewColOrbit = (CColOrbit<T> *)((char *)pNextRowColOrbit + nColCurr * colOrbLen);
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

	if (getUnforcedColOrbPntr()) {
        // Set unforced elements:
        for (auto row = firstUnforcedRow(); row <= nRow; row++) {
			const auto *pColOrbitIni = this->colOrbitIni(row);
			auto **ppUnforced = this->unforcedOrbits(row);
            for (int i = 1; i < this->rank(); i++) {
                pColOrbit = *(ppUnforced + i);
                while (pColOrbit) {
                    const size_t nColCurr = ((char *)pColOrbit - (char *)pColOrbitIni) / colOrbLen;
                    const auto lenFragm = pColOrbit->length();
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

template<class T>
void CEnumerator<T>::InitRowSolutions(const CEnumerator<T> *pMaster)
{
	const auto nRow = pMaster? pMaster->currentRowNumb() + 1 : 0;
	auto pSolutions = new CRowSolution<T>[rowNumb() - nRow];
	for (auto i = rowNumb(); i-- > nRow;)
		m_pRow[i] = pSolutions + i - nRow;

	if (nRow) 
		memcpy(rowStuffPntr(), pMaster->rowStuffPntr(), nRow * sizeof(*rowStuffPntr()));
}

template<class T>
size_t CEnumerator<T>::getDirectory(char *dirName, size_t lenBuffer) const
{
	const auto rowNumb = getInSys()->rowNumb();
#if 0
	const char *pDirNamePrefix = "V =";
	size_t lenPrefix = strlen(pDirNamePrefix);
	const size_t lenName = 4;
	const size_t lenDirName = lenPrefix + lenName + 1;
	if (lenBuffer <= lenDirName)
		return 0;

	// Define the length of text representation of the rowNumb in decimal format
	size_t len = lenName;
	auto nRow = rowNumb;
	while (nRow /= 10)
		len--;

	memcpy(dirName, pDirNamePrefix, lenPrefix);
	if (len > 0) {
		memset(dirName + lenPrefix, ' ', len);
		lenPrefix += len;
	}

	sprintf_s(dirName + lenPrefix, lenBuffer - lenPrefix, "%d", rowNumb);
#else
	const size_t lenDirName = sprintf_s(dirName, lenBuffer, "V =%4d", rowNumb);
#endif
	int retVal = 0;
	if (!fileExists(dirName, false))
		retVal = MKDIR(dirName);

	if (retVal)
		return 0;	// directory could not be used

	dirName[lenDirName] = '\\';
	return lenDirName + 1;
}

static bool getNextLineForComparison(FILE *file, char *buffer, int lenBuffer)
{
	bool outBlock = false;
	while (true) {
		if (!fgets(buffer, lenBuffer, file))
			return false;

		if (outBlock) {
			// Block was not open
			if (strstr(buffer, END_OUT_BLOCK))
				outBlock = false; // mark it as closed
		}
		else {
			// Block was not open
			outBlock = strstr(buffer, BEG_OUT_BLOCK) != NULL;
			if (!outBlock)
				return true;
		}
	}
}

template<class T>
bool CEnumerator<T>::compareResults(char *fileName, size_t lenFileName, bool *pBetterResults) const
{
	FOPEN(file, fileName, "r");
	if (!file)
		return false;

	// Create the name of the file with the previous results
	static size_t lenSuffix = strlen(FILE_NAME("")) + 1;
	strcpy_s(fileName + lenFileName, lenSuffix, FILE_NAME(""));

	FOPEN(filePrev, fileName, "r");
	if (!filePrev) {
		FCLOSE(file);
		return false;
	}

	bool retVal = false;			// the results are not the same
	const size_t lenBuf = 256;
	char buf[2][lenBuf];
	while (true) {
		if (!getNextLineForComparison(filePrev, buf[0], lenBuf))
			break;

		if (!getNextLineForComparison(file, buf[1], lenBuf))
			break;

		char *pntr0 = strstr(buf[0], CONSTRUCTED_IN);
		if (pntr0) {
			char *pntr1 = strstr(buf[1], CONSTRUCTED_IN);
			if (pntr1) {
				const size_t len = strlen(CONSTRUCTED_IN);
				*pBetterResults = CEnumInfo<T>::compareTime(pntr0 + len, pntr1 + len);
				retVal = true;		// the results are the same
			}

			break;
		}

		if (strcmp(buf[0], buf[1]))
			break;
	}

	FCLOSE(file);
	FCLOSE(filePrev);
	return retVal;
}
/*
template<class T>
bool CEnumerator<T>::printMatrix(const designRaram *pParam) const
{
	const uint outType = pParam->outType;
	return	outType & t_AllObject ||
		outType & t_Transitive && groupIsTransitive() ||
		outType & t_GroupOrderGT && groupOrder() > pParam->grpOrder ||
		outType & t_GroupOrderLT && groupOrder() < pParam->grpOrder ||
		outType & t_GroupOrderEQ && groupOrder() == pParam->grpOrder;
}
*/
#if USE_STRONG_CANONICITY_A
void CEnumerator<T>::checkUnusedSolutions(CRowSolution *pRowSolution)
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
		const auto *pCurrSol = pRowSol->currSolution();
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

#if CANON_ON_GPU
template<class T>
size_t CEnumerator<T>::copyColOrbitInfo(T nRow) const
{
	// Function copy existing information about the orbits of columns  
	// into new structures, which could be used on GPU
	const auto pColOrbInfoBeg = GPU_CanonChecker()->ColOrbitData(t_CPU);
	auto pColOrbInfo = pColOrbInfoBeg;
	const auto colOrbit = colOrbits();
	const auto colOrbitIni = colOrbitsIni();
	while (nRow--) {
		const auto pColOrbitIni = colOrbitIni[nRow];
		auto pColOrbit = colOrbit[nRow];
		if (!pColOrbit) {
			*pColOrbInfo++ = UINT64_MAX;
			continue;
		}

		// Keep pointer to set number of elements stored in following loop
		auto pColOrbInfoNum = pColOrbInfo++;
		*pColOrbInfo++ = (char *)pColOrbit - (char *)pColOrbitIni;
		while (pColOrbit) {
			*pColOrbInfo++ = pColOrbit->length();
			pColOrbit = pColOrbit->next();
			*pColOrbInfo++ = (char *)pColOrbit - (char *)pColOrbitIni;
		}

		*pColOrbInfoNum = ((--pColOrbInfo - pColOrbInfoNum) >> 1);
	}

	return pColOrbInfo - pColOrbInfoBeg;
}

#if TRACE_CUDA_FLAG
static int cntr;
#endif

#endif
