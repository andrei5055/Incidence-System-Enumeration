
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
	#define MKDIR(x)  mkdir(dirName, 0777)
#else
	"Please define corresponding method for making directory"
#endif

#include <mutex>  

#if PRINT_SOLUTIONS || PRINT_CURRENT_MATRIX
int ccc = 0;
bool startPrinting = START_PRINTING_AFTER <= 0;
#endif

#if USE_MUTEX
std::mutex out_mutex;
#endif

template class CEnumerator<TDATA_TYPES>;

FClass2(CEnumerator, bool)::fileExists(const char *path, bool file) const
{
	struct stat info;
	return (stat(path, &info) == 0) && (file || info.st_mode & S_IFDIR);
}

FClass2(CEnumerator, RowSolutionPntr)::FindRowSolution(S *pPartNumb)
{
	// It's OK to use permStorage() and not permStorage(nPart) here,
	// because the values permStorage(nPart) are the same for all nPart
	this->setUseCanonGroup(USE_CANON_GROUP && !permStorage()->isEmpty());

	const auto numParts = this->numParts();
	auto * const pSolutionWereConstructed = getSolutionsWereConstructed(numParts, currentRowNumb());
	auto i = *pPartNumb;
	if (i) {
		// Only for combined designs
		// Check if the solutions for current row were already constructed before for previous parts
		while (i-- && !pSolutionWereConstructed[i]);
		i++;
	}

	S nVar = ELEMENT_MAX;
	const auto firstPart = i;
	if (prepareToFindRowSolution()) {
		// Find row solution for all parts of the design
		while (true) {
			if (pSolutionWereConstructed && i >= firstPart)
				pSolutionWereConstructed[i] = 0;

			nVar = MakeSystem(i);
			if (nVar == ELEMENT_MAX)
				break;    // There is no system of equations which could have valid solutions for next row

			if (i < firstPart) {
				i++;
				continue; // The solutions up to firstPart where already constructed
			}

			// NOTE: For incidence system nVar could be 0 for last row, when on current row none of the orbits was splitted into two parts
			setPrintResultNumVar(nVar);
			RowSolutionPntr pRowSolution = FindSolution(nVar, i, m_lastRightPartIndex[i]);
			if (!pRowSolution || !pRowSolution->numSolutions())
				break;

			if (!checkSolutions(pRowSolution, i, m_lastRightPartIndex[i]))
				break;

			if (pSolutionWereConstructed)
				pSolutionWereConstructed[i] = 1;

			if (numParts > 1) {
				if (i)
					pRowSolution->setSolutionIndex(0);
				else
					pRowSolution->saveSolutionIndex();
			}

			if (++i >= numParts)
				break;

			if (i == 1) {
				// Save index to the first solution to be tested for part 0 
				setFirstPartSolutionIndex(pRowSolution->solutionIndex());
			}
		}
	}

	OUTPUT_SOLUTION(this->rowStuff(currentRowNumb()), outFile(), currentRowNumb(), false, firstPart, i);
	if (i < numParts) {
		*pPartNumb = i;
#if PRINT_SOLUTIONS
		if (MAKE_OUTPUT()) {
			char buff[128], * pBuff = buff;
			pBuff += SNPRINTF(buff, sizeof(buff), "\nNo solutions for ");
			if (numParts > 1)
				pBuff += SNPRINTF(pBuff, sizeof(buff) - (pBuff - buff), "the part %d of ", i);
			SNPRINTF(pBuff, sizeof(buff) - (pBuff - buff), "the row %d. System of equations %s\n", currentRowNumb(),
				nVar == ELEMENT_MAX ? "cannot be constructed" : "has no valid solutions");
			outString(buff, outFile());
		}
#endif
		return NULL;
	}

	// Always return the pointer to the solutions of part #0
	return this->rowStuff(currentRowNumb());
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

template<typename T, typename S>
void threadEnumerate(Class2(CThreadEnumerator) *threadEnum, designParam *param, const EnumeratorPntr pMaster)
{
	threadEnum->EnumerateBIBD(param, pMaster);
}

FClass2(CEnumerator, int)::threadWaitingLoop(int thrIdx, t_threadCode code, Class2(CThreadEnumerator) *threadEnum, size_t nThread) const
{
	// Try to find the index of not-running thread
	const auto flag = code == t_threadNotUsed;
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

		auto *pEnum = threadEnum + thrIdx;
#if WRITE_MULTITHREAD_LOG 
		fprintf(file, "==>thrIdx = %2d  CODE = %d  nThread= %d\n", thrIdx, pEnum->code(), nThread);
		fclose(file);
#endif
		switch (pEnum->code()) {
			case t_threadFinished:
				LAUNCH_CANONICITY_TESTING(pEnum, this);
				this->enumInfo()->updateGroupInfo(pEnum->enumInfo());
				thread_message(thrIdx, "finished", code);
				if (flag) {
					// Changing code, otherwise we in next call CEnumInfo<T>::reportProgress(...)
					// t_canonical/t_totalConstr matrices will be counted one more time
					pEnum->setCode(t_threadNotUsed); 
				}

				this->enumInfo()->reportProgress(pEnum);
				if (flag)
					continue;

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

FClass2(CEnumerator, bool)::Enumerate(designParam* pParam, bool writeFile, EnumInfoPntr pEnumInfo, const EnumeratorPntr pMaster, t_threadCode* pTreadCode)
{
	setDesignParams(pParam);
	const auto* pInpMaster = pMaster;
#if !CONSTR_ON_GPU
	static std::mutex mtx;
	char buff[256], jobTitle[256];
	const auto lenBuffer = countof(buff);
	const auto threadNumb = pParam->threadNumb;
	// We will not launch separate thread, when threadNumb is equal to 1
	const int mt_level = threadNumb >= 1? pParam->mt_level : INT_MAX;

	size_t lenName = 0;
	bool knownResults = false;
	if (writeFile) {
		// We will be here only for the master
		pParam->firstMatr = true;
		if (!makeJobTitle(pParam, jobTitle, countof(jobTitle), "\n"))
			return false;

		// Construct the file name of the file with the enumeration results
		if (!getMasterFileName(buff, lenBuffer, &lenName))
			return false;

		// The results are known, if the file with the enumeration results exists
		knownResults = fileExists(buff);
		if (knownResults)
			strcpy_s(buff + lenName, lenBuffer - lenName, FILE_NAME(CURRENT_RESULTS));
		else
			lenName = 0;

		// Create a new file for output of the enumeration results
		const auto newFile = this->createNewFile(buff);
		const auto seekFile = !newFile && SeekLogFile() ? long(designParams()->rewindLen) : 0;
		FOPEN(file, buff, newFile ? "w" : (seekFile ? "r+t" : "a"));
		if (seekFile)  // Previously written end of the log file needs to be removed
			fseek(file, -seekFile, SEEK_END);

		this->setOutFile(file);
		outString(jobTitle, this->outFile());
	} else
		this->setOutFile(NULL);

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
	ThreadEnumeratorPntr pThreadEnum = NULL;
	int thrIdx = 0;
	#if USE_POOL
		asio::io_service *pIoService = NULL;
		thread_group *pThreadpool = NULL;
	#endif
#endif
	const auto pMatrix = this->matrix();

	const auto firstNonfixedRow = firtstNonfixedRowNumber();
	// Allocate memory for the orbits of two consecutive rows
	T lenStab;
	RowSolutionPntr pRowSolution;
	RowSolutionPntr pNextRowSolution = NULL;
	InitRowSolutions(pMaster);
	const auto nRows = rowNumb();
	const bool threadFlag = pMaster != NULL;
	if (threadFlag) {
		lenStab = pMaster->stabiliserLengthExt();
		this->setCurrentRowNumb(nRow = pMaster->currentRowNumb());
		if (pMaster->getGroupOnParts())
			updateCanonicityChecker(rowNumb(), colNumb());

		pRowSolution = rowStuff(nRow);
		InitGroupOderStorage(pMaster->getGroupOnParts());
		this->setOutFile(pMaster->outFile());
		setX0_3(pMaster->getX0_3());
		const auto firstUnforced = pMaster->firstUnforcedRow();
		if (firstUnforced > 0) {
			setFirstUnforcedRow(firstUnforced);

			// Unforced lambda's are stored one by one for all parts and rows
			const auto length = numParts() * (nRows - firstUnforced) * sizeof(*forcibleLambdaPntr());
			const auto from = firstUnforced * numParts();
			memcpy(forcibleLambdaPntr() + from, pMaster->forcibleLambdaPntr() + from, length);
		}
	} else {
		lenStab = nRow = firstNonfixedRow - 2;
		CreateForcedRows();
		pRowSolution = setFirstRowSolutions();
		this->setEnumInfo(pEnumInfo);
		pEnumInfo->startClock();

#if USE_THREADS_ENUM 
		if (pParam->threadNumb)
			pThreadEnum = new Class2(CThreadEnumerator)[pParam->threadNumb];
		CMatrixData<TDATA_TYPES>::ResetCounter();
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
    
	this->setStabiliserLengthExt(lenStab);
#if CANON_ON_GPU
	this->matrix()->m_nStabExtern = lenStab;
#endif
	if (!threadFlag)
		prepareToTestExtraFeatures();

	// For multi-threaded version we need to test only one top level solution
	const bool multiPartDesign = numParts() > 1;
	const auto outInfo = t_saveRowToChange + t_saveRowPermutations;
	const auto nRowEnd = nRow ? nRow + 1 : 0;

	const auto use_master_sol = designParams()->use_master_sol;
	this->initiateColOrbits(nRows, nRow, pMatrix->partsInfo(), this->IS_enumerator(), use_master_sol, pMaster);

	// Construct nontrivial group, acting on parts (groups of blocks)
	// As of today (09/29/2021), it should work only for CombBIBD's with at least two equal lambda's
	auto* pGroupOnParts = pMaster ? pMaster->getGroupOnParts() : makeGroupOnParts(this);
	setGroupOnParts(pGroupOnParts);
	MatrixDataPntr pSpareMatrix = pGroupOnParts? CreateSpareMatrix(pMaster) : NULL;

	T level(0), nPart(0);
	TestCanonParams<T,S> canonParam = {this, matrix(), numParts(), &nPart, &level, pGroupOnParts, pSpareMatrix};

	CreateAuxiliaryStructures(pMaster);

	const auto procFlag = designParams()->find_master_design || designParams()->find_all_2_decomp;
	// minimal index of the part, which will be changed on current row
	auto *firstPartIdx = new T[nRows];
	memset(firstPartIdx, 0, nRows * sizeof(*firstPartIdx));
	bool canonMatrix = false;

	while (pRowSolution) {
		const bool useCanonGroup = USE_CANON_GROUP && nRow > 0;
		bool checkNextPart = false;
		auto iFirstPartIdx = numParts();
#if USE_THREADS_ENUM
		const bool usingThreads = pThreadEnum && nRow == mt_level;
		if (usingThreads) {
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
						thread t1(threadEnumerate<T, S>, pThreadEnum + thrIdx, pParam, this);
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
#else
			const bool usingThreads = false;
#endif
			REPORT_PROGRESS(pEnumInfo, t_reportCriteria::t_reportByTime);
			OUTPUT_SOLUTION(pRowSolution, outFile(), nRow, true, 0, numParts());
			MakeRow(pRowSolution, nRow == firstNonfixedRow, firstPartIdx[nRow]);

			canonMatrix = true;
			if (++nRow == nRows) {
				if (multiPartDesign)
					pRowSolution->restoreSolutionIndex();  // Andrei: Perhaps for CombBIBD with nPart > 2 we need to do it a more sophisticated way

				pEnumInfo->incrConstrTotal();
				bool flag = true;

				if (!TestCanonicityOnGPU()) {
					EXIT(-1);
					canonMatrix = this->TestCanonicity(nRow, &canonParam, outInfo);
					if (canonMatrix) {
						if (procFlag)
							ConstructedDesignProcessing();

						//	DEBUGGING: How Construct Aut(D): int ddd = canonChecker()->constructGroup();
						int matrFlags = 0;
						if (TestFeatures(pEnumInfo, pMatrix, &matrFlags, pMaster)) {
							if (noReplicatedBlocks() && pEnumInfo->constructedAllNoReplBlockMatrix()) {
								pEnumInfo->setNoReplBlockFlag(false);
								level = getInSys()->GetK();
								flag = false;
							}
							else {
								pEnumInfo->updateConstrCounters(matrFlags, this);
#if !CONSTR_ON_GPU
								if (this->printMatrix(pParam)) {
									mtx.lock();
									outBlockTitle();
#if USE_THREADS
									pMatrix->printOut(this->outFile(), nRow, 0, this);
#else
									pMatrix->printOut(this->outFile(), nRow, pEnumInfo->numMatrOfType(t_design_type::t_canonical), this);
#endif
									mtx.unlock();
								}
#endif
#if PRINT_SOLUTIONS
								ccc = 0;
#endif
								if (!this->rowMaster())  // We are not in the slave thread
									REPORT_PROGRESS(pEnumInfo, t_reportCriteria::t_matrConstructed);
							}
						}
					}
					else {
						OUTPUT_MATRIX(pMatrix, outFile(), nRow, pEnumInfo, false);
						if (!multiPartDesign)  // NOTE: as of today (08/27/2021) "level" is set correctly only when when nPart = 0
							flag = false;
					}
				}

				if (!flag) {
					if (outInfo & t_saveRowToChange) {
						while (--nRow > level)
							this->reset(nRow);
					}
					else
						--nRow;

					pRowSolution = rowStuff(nRow);
					this->setCurrentRowNumb(nRow);
				}
				else {
					checkNextPart = multiPartDesign;
					nRow--;
					if (canonMatrix && multiPartDesign) {
						pRowSolution = rowStuff(nRow-1);
						while (--iFirstPartIdx) {
							auto pPartRowSolution = pRowSolution + iFirstPartIdx;
							this->resetUnforcedColOrb(iFirstPartIdx);
							if (pPartRowSolution->allSolutionChecked())
								pPartRowSolution->setSolutionIndex(0);
							else
								break;
						}
						if (iFirstPartIdx) {
							this->setCurrentRowNumb(--nRow);
							firstPartIdx[nRow] = iFirstPartIdx;
							continue;
						}
						this->resetUnforcedColOrb(0);
					}

					pRowSolution = NULL;
				}
			}
			else {
				this->setCurrentRowNumb(nRow);
				for (auto i = numParts(); i-- > firstPartIdx[nRow-1];) {
					this->setColOrbitCurr(m_pFirstColOrb[i], i);
					this->setCurrUnforcedOrbPtr(nRow, i);
				}

				canonMatrix = !USE_CANON_GROUP || 
					          this->TestCanonicity(nRow, &canonParam, t_saveNothing, pRowSolution);
				OUTPUT_MATRIX(pMatrix, outFile(), nRow, pEnumInfo, canonMatrix);

				if (canonMatrix) {
					if (pInpMaster) {
						copyInfoFromMaster(pInpMaster);
#if WAIT_THREADS
						pInpMaster = NULL;
						*pTreadCode = t_threadRunning;
#endif
					}

					if (!useCanonGroup)
						this->setGroupOrder(1);

					setPrintResultRowNumber(nRow);
					pRowSolution = FindRowSolution(firstPartIdx + nRow - 1);
#if USE_THREADS && !WAIT_THREADS
					if (pInpMaster && pRowSolution) {
						pInpMaster = NULL;
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
					if (multiPartDesign) {
						// When solution for the first part is not canonical, we don't have
						// to check all combinations of solutions for remaining parts
						const auto changeFirstPart = !*canonParam.pPartNumb;
						while (--iFirstPartIdx) {
							auto pPartRowSolution = pRowSolution + iFirstPartIdx;
							this->resetUnforcedColOrb(iFirstPartIdx, nRow);
							if (changeFirstPart || pPartRowSolution->allSolutionChecked())
								pPartRowSolution->setSolutionIndex(0);
							else
								break;
						}

						if (iFirstPartIdx) {
							// Enumeration of combined designs AND not all solutions for i-th part were tested
							this->setCurrentRowNumb(--nRow);
							firstPartIdx[nRow] = iFirstPartIdx;
							continue;
						}

						// We will change the first portion of CombBIB. Therefore,
						// everything which was enforced by this portion should be reset.
						this->resetUnforcedColOrb(0, nRow);   
					}

					checkNextPart = multiPartDesign;
					pRowSolution = NULL;
				}
			}
#if USE_THREADS_ENUM
		}
#endif
		if (multiPartDesign) {
			// We enumerate the multi-part designs
			if (usingThreads) {
				// We just finished with all solutions tested by the threads.
				rowStuff(nRow)->restoreSolutionIndex();
			} else
			if (!pRowSolution) {
				if (iFirstPartIdx) {
					// We can reach this point by three different paths:
					// 1. We are not using threads AND
					//   (a) matrix is NOT canonical OR
					//   (b) matrix was canonical, but solution for one of the parts was not found
					// 2. We are using threads AND we are done with next rows.
					// When we do have (1.b), but firstPartIdx[nRow] == 0, we should proceed as in case (1.a)
					--nRow;
					auto lastPartIdx = usingThreads? numParts()-1 : firstPartIdx[nRow];
					bool check_all_solution = !canonMatrix;
					if (canonMatrix) {
						// When there is no solution for some part, the indices
						// for all remaining parts are equal to 0
						pRowSolution = rowStuff(nRow);
						while (lastPartIdx && (pRowSolution + lastPartIdx)->allSolutionChecked()) {
							this->resetUnforcedColOrb(lastPartIdx, nRow + 1);
							(pRowSolution + lastPartIdx--)->setSolutionIndex(0);
						}

						this->resetUnforcedColOrb(lastPartIdx, nRow + 1);
						if (lastPartIdx) {
							firstPartIdx[nRow] = lastPartIdx;
							this->setCurrentRowNumb(nRow);
							if (usingThreads) {
								// We need to restart the whole process on the threads for next row
								rowStuff(nRow + 1)->restoreSolutionIndex();
							}
							continue;
						}
					}

					iFirstPartIdx = numParts() - 1;
					while (true) {
						pRowSolution = rowStuff(nRow, iFirstPartIdx);
						if (check_all_solution && !pRowSolution->allSolutionChecked()) {
							firstPartIdx[nRow] = iFirstPartIdx;
							break;
						}

						pRowSolution->setSolutionIndex(0);
						if (!--iFirstPartIdx)
							break;

						if (iFirstPartIdx == lastPartIdx)
							check_all_solution = true;
					}

					this->setCurrentRowNumb(nRow);
					pRowSolution = rowStuff(nRow);

					for (auto i = numParts(); i-- > iFirstPartIdx;)
						this->resetUnforcedColOrb(i, nRow + 1);

					if (iFirstPartIdx)
						continue;
				}
				else {
					pRowSolution = rowStuff(--nRow);
					setCurrentRowNumb(nRow);
					for (auto i = numParts(); --i;)
						(pRowSolution+i)->setSolutionIndex(0);
				}
			}

			if (nRow < nRowEnd)
				break;

			bool resetSolutions = true;
			firstPartIdx[nRow] = 0;  // When we are here the first fragment of CombBIB will be changed
			while (!(pNextRowSolution = pRowSolution) || !(pRowSolution = pRowSolution->NextSolution(useCanonGroup))) {
				this->reset(nRow, resetSolutions);
				const auto nRowNext = nRow--;
				if (pNextRowSolution)
					pNextRowSolution->restoreSolutionIndex();

				pRowSolution = rowStuff(nRow);
				auto j = numParts();
				const auto* ppOrb = colOrbitPntr() + nRow * pMatrix->colNumb();
				const auto partsInfo = pMatrix->partsInfo();
				while (--j) {
					auto pPartRowSolution = pRowSolution + j;
					this->resetUnforcedColOrb(j, nRowNext);    // New code
					this->setCurrUnforcedOrbPtr(nRow, j);
					if (!pPartRowSolution->allSolutionChecked())
						break;

					setColOrbitCurr(*(ppOrb + partsInfo->getShift(j)), j);
					pPartRowSolution->setSolutionIndex(0);
				}


				if (firstPartIdx[nRow] = j) {
					this->setCurrentRowNumb(nRow);
					firstPartIdx[nRow] = 1;
					break;
				}

				if (nRow < nRowEnd) {
					pRowSolution = NULL;
					break;
				}

				rowStuff(nRowNext)->resetSolution();
				ppOrb = colOrbitPntr() + nRowNext * pMatrix->colNumb();
				setColOrbitCurr(*ppOrb, 0);
				this->resetUnforcedColOrb(0);   // New code 
				this->resetFirstUnforcedRow();
				this->setCurrentRowNumb(nRow);
			}
		}
		else {
			while (!pRowSolution || !(pRowSolution = pRowSolution->NextSolution(useCanonGroup))) {
				this->reset(nRow);
				if (nRow-- <= nRowEnd)
					break;

				setCurrentRowNumb(nRow);
				pRowSolution = rowStuff(nRow);
			}
		}
	} 

	delete[] firstPartIdx;

#if USE_THREADS_ENUM && !WAIT_THREADS
	if (pThreadEnum)
		threadWaitingLoop(thrIdx, t_threadNotUsed, pThreadEnum, pParam->threadNumb);
#endif

    this->closeColOrbits(use_master_sol);

	if (!threadFlag || !USE_THREADS_ENUM) {

#if USE_THREADS_ENUM
		delete[] pThreadEnum;
#else
#ifdef USE_CUDA
		// This method is called after thread is ended, When they are used
		LAUNCH_CANONICITY_TESTING(enumInfo(), this);
#endif
#endif

#if !CONSTR_ON_GPU
		// We are not in the slave thread
		if (!pParam->firstMatr)
			outString("\n" END_OUT_BLOCK "Constructed Matrices " BEG_OUT_BLOCK "\n", this->outFile());

		beforeEnumInfoOutput();
		pEnumInfo->outEnumInfo(this->outFilePntr(), lenName == 0);

		closeFile();
		t_resType resType;
		if (lenName) {
			// Compare current results with previously obtained
			bool betterResults = true;
			const char *currentFile = FILE_NAME(CURRENT_RESULTS);
			strcpy_s(buff + lenName, countof(buff) - lenName, currentFile);
			// TO DO: For Semi-Symmetric graphs more complicated comparison function should be implemented
			if (pParam->objType != t_objectType::t_SemiSymmetricGraph && compareResults(buff, lenName, &betterResults)) {
				// Create the file name with the current results 
				strcpy_s(jobTitle, countof(jobTitle), buff);
				strcpy_s(jobTitle + lenName, countof(jobTitle) - lenName, currentFile);

				if (betterResults) {
					remove(buff);			// Remove file with previous results
					rename(jobTitle, buff);	// Rename file
					resType = t_resType::t_resBetter;
				}
				else {
					if (pParam->firstMatr)
						remove(jobTitle);	// Deleting new file only when it does not contain matrices and 

					resType = t_resType::t_resWorse;
				}
			}
			else
				resType = t_resType::t_resInconsistent; // results are not the same as before
		}
		else {
			if (pParam->objType != t_objectType::t_SemiSymmetricGraph) {
				if (getMasterFileName(buff, lenBuffer, &lenName))
				    compareResults(buff, lenName);
			}
			resType = t_resType::t_resNew;
		}

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
			*pTreadCode = pInpMaster ? t_threadLaunchFail : t_threadFinished;
	}

	return true;
} 

FClass2(CEnumerator, void)::reset(T nRow, bool resetSolutions) {
	const auto *pMatrix = this->matrix();
	const auto *ppOrb = colOrbitPntr() + pMatrix->colNumb() * nRow;

	auto j = numParts();
	if (j > 1) {
		const auto partsInfo = pMatrix->partsInfo();
		while (--j) {
			setColOrbitCurr(*(ppOrb + partsInfo->getShift(j)), j);
			this->resetUnforcedColOrb(j);  // New Code
		}
	}
	else {
		rowStuff(nRow)->resetSolution();
		setColOrbitCurr(*ppOrb, 0);
		this->resetUnforcedColOrb(0);       // Regular code
	}

	this->resetFirstUnforcedRow();
}

FClass2(CEnumerator, void)::MakeRow(RowSolutionPntr pRowSolution, bool flag, S iFirstPartIdx) {
	flag &= !iFirstPartIdx;   // X0_3 condition should be checked only on first part
	                          // of CombBIBD or on regular incidence system
	// Loop over all portions of the solution
	const auto nRow = currentRowNumb();
	auto* const pSolutionWereConstructed = getSolutionsWereConstructed(numParts(), nRow + 1);
	for (auto i = 0; i < numParts(); i++) {
		auto pPartRowSolution = pRowSolution + i;
		// We need to get lastRightPartIndex here and use later because 
		// for multi-thread configuration it could be changed by master
		if (i) {
			// When we are in that function, the solutions for the first part was just changed
			// It means that we need to check all combinations of solutions for remaining parts
			m_lastRightPartIndex[i] = pPartRowSolution->numSolutions() - 1;
		} else
			m_lastRightPartIndex[i] = pPartRowSolution->solutionIndex();

		if (i < iFirstPartIdx)
			continue;

		const auto pCurrSolution = pPartRowSolution->currSolution();
		m_pFirstColOrb[i] = CMatrixCanonChecker::MakeRow(nRow, pCurrSolution, true, i);

		if (flag) {
			flag = false;
			setX0_3(*pCurrSolution);
		}

		if (pSolutionWereConstructed)           // All previously constructed solutions for current portion (if any) are not valid
			pSolutionWereConstructed[i] = 0;	// anymore because we just changed corresponding fragments of previous row.
	}
}

FClass2(CEnumerator, void)::InitRowSolutions(const EnumeratorPntr pMaster) {
	const auto use_master_sol = designParams()->use_master_sol;
	const auto nRow = pMaster? pMaster->currentRowNumb() + use_master_sol : 0;
	const auto pMatrix = pMaster? pMaster->matrix() : this->matrix();
	const auto nParts = pMatrix->numParts();
	auto i = rowNumb();
	const auto pSolutions = new Class2(CRowSolution)[nParts * (i - nRow)];
	while (i-- > nRow)
		m_pRow[i] = pSolutions + nParts * (i - nRow);

	if (pMaster) {
		if (!use_master_sol) {
			const auto* pRowSolution = pMaster->rowStuff(nRow);
			for (auto i = 0; i < numParts(); i++)
				*(pSolutions + i) = *(pRowSolution + i);
		}

		// Just in case, copying pointers to the solutions from master
		memcpy(rowStuffPntr(), pMaster->rowStuffPntr(), nRow * sizeof(*rowStuffPntr()));
	}
}

FClass2(CEnumerator, size_t)::getDirectory(char *dirName, size_t lenBuffer, bool rowNeeded) const
{
	const auto pParam = designParams();
	lenBuffer--;		// Reserving 1 byte for last '/'

	auto len = SNPRINTF(dirName, lenBuffer, "%s", pParam->workingDir.c_str());
	if (fileExists(dirName, false) ? 0 : MKDIR(dirName))
		return 0;		// Directory could not be used

	const auto* pDirName = this->getTopLevelDirName();
	if (pDirName) {
		len += SNPRINTF(dirName + len, lenBuffer - len, "%s/", pDirName);
		if (fileExists(dirName, false) ? 0 : MKDIR(dirName))
			return 0;	// Directory could not be used
	}

	if (rowNeeded) {
		auto rowNumb = getInSys()->rowNumbExt();
		if (pParam->objType == t_objectType::t_SemiSymmetricGraph)
			rowNumb *= pParam->r / pParam->k;

		len += SNPRINTF(dirName + len, lenBuffer - len, "V =%4d/", rowNumb);
		if (fileExists(dirName, false) ? 0 : MKDIR(dirName))
			return 0;	// Directory could not be used
	}

	return len;
}

static bool getNextLineForComparison(FILE *file, char *buffer, int lenBuffer, char *tmpBuff)
{
	// characters to skip at the beginning of each line
	bool skipTable[256] = {};
	skipTable[' '] = skipTable['\n'] = skipTable['\t'] = true;

	bool outBlock = false;
	while (true) {
		if (!fgets(tmpBuff, lenBuffer, file))
			return false;

		if (outBlock) {
			// Block was not open
			if (strstr(tmpBuff, END_OUT_BLOCK))
				outBlock = false; // mark it as closed
		}
		else {
			// Block was not open
			outBlock = strstr(tmpBuff, BEG_OUT_BLOCK) != NULL;
			if (!outBlock) {
				if (strstr(tmpBuff, OF_THEM))
					continue;

				// Remove leading spaces:
				auto* pCurrentSymb = tmpBuff;
				while (skipTable[*pCurrentSymb++]);
				if (!*(--pCurrentSymb)) // Is it the end of line?
					continue;			// Skip an "empty" line

				strncpy_s(buffer, lenBuffer, pCurrentSymb, lenBuffer);
				return true;
			}
		}
	}
}

FClass2(CEnumerator, bool)::getMasterFileName(char *buffer, size_t lenBuffer, size_t *pLenName) const
{
	// Construct the file name of the file with the enumeration results
	if (!makeFileName(buffer, lenBuffer, ""))
		return false;

	*pLenName = strlen(buffer);

	// Create file name for the output of final enumeration results
	strcpy_s(buffer + *pLenName, lenBuffer - *pLenName, FILE_NAME(""));
	return true;
}

FClass2(CEnumerator, bool)::cmpProcedure(FILE* file[2], bool *pBetterResults)
{
	const size_t lenBuf = 256;
	const size_t len = strlen(CONSTRUCTED_IN);
	char buf[3][lenBuf] = {}, * pntr[2] = { NULL, NULL };
	while (true) {
		bool eof = false;
		for (int i = 0; i < 2; i++) {
			if (file[i]) {
				if (!getNextLineForComparison(file[i], buf[i], lenBuf, buf[2])) {
					pntr[i] = 0;
					if (!i)
						continue;

					return true;
				}

				pntr[i] = strstr(buf[i], CONSTRUCTED_IN);
			}
		}

		if (pntr[0]) {
			if (!file[1] || pntr[1]) {
				if (pBetterResults)
					*pBetterResults = Class2(CEnumInfo)::compareTime(pntr[0] + len, pntr[1] + len);

				if ((!pBetterResults || *pBetterResults) && getNextLineForComparison(file[0], buf[1], lenBuf, buf[2])) {
					char* pInfo[] = { pntr[0] + len, buf[0], buf[1] };
					UpdateEnumerationDB(pInfo, 3);
				}

				return true;		// the results are the same
			}

			break;
		}

		if (file[1] && strcmp(buf[0], buf[1])) {
			if (buf[0][0] == '_' && buf[0][0] == buf[1][0])
				continue;  // skip decorating line

			// Try to compare each word separately
			// Trim leading "new lines"
			size_t len[] = { strlen(buf[0]), strlen(buf[1]) };
			for (int i = 0; i < 2; i++) {
				while (*(buf[i] + len[i] - 1) == '\n')
					len[i]--;
			}

			char* pBeg[2] = {}, * pEnd[] = { buf[0] - 1, buf[1] - 1 };
			int flag = 0;
			while (!flag) {
				const bool cond = buf[0] + len[0] == pEnd[0];
				if (cond != (buf[1] + len[1] == pEnd[1])) {
					flag = 1;
					break;
				}
				else
				if (cond)
					break;

				for (int i = 0; i < 2; i++) {
					pBeg[i] = pEnd[i];
					while (*(++pBeg[i]) == ' ');
					char* p = pBeg[i];
					while (*++p != ' ' && *p != '\n' && *p != '\0');
					*(pEnd[i] = p) = '\0';
				}
				flag = strcmp(pBeg[0], pBeg[1]);
			}

			if (flag)
				break;
		}
	}

	return false;
}

FClass2(CEnumerator, bool)::compareResults(char *fileName, size_t lenFileName, bool *pBetterResults) {
	FOPEN(file, fileName, "r");
	if (!file)
		return false;

	FILE* ppFile[] = { file, NULL };

	if (pBetterResults) {
		// Create the name of the file with the previous results
		static auto lenSuffix = strlen(FILE_NAME("")) + 1;
		strcpy_s(fileName + lenFileName, lenSuffix, FILE_NAME(""));
		FOPEN(filePrev, fileName, "r");
		ppFile[1] = filePrev;
	}

	const bool retVal = cmpProcedure(ppFile, pBetterResults);
	FCLOSE(file);
	FCLOSE(ppFile[1]);

	return ppFile[1]? retVal : false;
}

static void outKeyInfo(const char* key, char **pInfo, FILE* file, const char *pComment = NULL)
{
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	localtime_s(&tstruct, &now);
	strftime(buf, sizeof(buf), "%b %d, %Y", &tstruct);
	fprintf(file, "%s  %12s   %9s  %15s  %13s", key, pInfo[1], pInfo[2], pInfo[0], buf);
	if (pComment)
		fprintf(file, "%s", pComment);
	else
		fprintf(file, "\n");
}

FClass2(CEnumerator, void)::outputTitle(FILE *file) const {
	char format[64];
	SPRINTF(format, "%s%%9s:  %%9s:       %%9s: %%9s:      %%9s:\n", this->getObjNameFormat());
	fprintf(file, format, this->getTopLevelDirName(), "Total #", "Simple #", "Run Time", "Date", "Comments");
}

FClass2(CEnumerator, void)::UpdateEnumerationDB(char **pInfo, int len)
{
	for (int i = 0; i < len; i++) {
		// Eliminate all first spaces
		auto* pntr = pInfo[i];
		while (*pntr == ' ')
			++pntr;

		int j = 0;
		auto* pEnd = strstr(pInfo[i] = pntr, "\n");
		if (pEnd)
			*pEnd = '\0';

		while (*pntr && *pntr != ' ')
			pntr++;

		*pntr = '\0';
	}

	char enumerationDB[256];
	const auto length = getDirectory(enumerationDB, countof(enumerationDB), false);
	sprintf_s(enumerationDB + length, countof(enumerationDB) - length, "EnumerationDB.txt");

	FILE *dbFile = NULL;
	int i = -1;
	while (!dbFile && ++i < 2) {
		FOPEN(file, enumerationDB, i? "w" : "r");
		dbFile = file;
	}

	if (!dbFile)
		return;

	char key[32], adjustedKey[64];// , keyCmp[64];
	this->getEnumerationObjectKey(key, countof(key));
	// the lexicografical order of key's could be adjusted for some type of designs using:
	const char *pAdjKey = this->getEnumerationObjectKeyA(adjustedKey, countof(adjustedKey)/2);

	if (i) {
		this->outputTitle(dbFile);
		outKeyInfo(key, pInfo, dbFile);
		fclose(dbFile);
		return;
	}

	int resCmp = -1;
	char tmpFile[256], buffer[256];
	// EnumerationDB file exists
	// We need to find a record which corresponds to the key
	sprintf_s(tmpFile, "%s_tmp", enumerationDB);
	FOPEN(f, tmpFile, "w");
	bool compareFlag = true;
	const auto lenKey = strlen(key);
	const auto lenKeyAdj = pAdjKey? strlen(pAdjKey) : 0;
	bool firstLine = true;   // We don't need to compare key with the first line
	while (fgets(buffer, countof(buffer), dbFile)) {
		if (buffer[0] != ';') {
			// Not a comment line
			if (compareFlag && !firstLine) {
				if (buffer[0] && !strchr("_-", buffer[0])) {
					resCmp = pAdjKey? compareEnumerationDB_record(buffer) : strncmp(buffer, key, lenKey);
				}
				else
					resCmp = 1;   // EOF found

				if (resCmp >= 0) {
					const char* pComment = !resCmp? strstr(buffer, " >> ") : NULL;
					outKeyInfo(key, pInfo, f, pComment);
					compareFlag = false;
					if (!resCmp)
						continue;
				}
			}
		}

		firstLine = false;
		if (fputs(buffer, f) < 0) {
			FCLOSE(dbFile);
			FCLOSE(f);
			return; // Something wrong
		}
	}

	if (compareFlag)
		outKeyInfo(key, pInfo, f);

	FCLOSE(dbFile);
	FCLOSE(f);

	std::string error;
	if (remove(enumerationDB) != 0) {
		error = " Cannot remove file '";
	}
	else
	if (rename(tmpFile, enumerationDB) != 0) {
		error = " Cannot rename file '";
		error += tmpFile;
		error += "' to '";
	}
	else
		return;

	error += enumerationDB;
	error += "'";
	perror(error.c_str());
}

FClass2(CEnumerator, void)::outBlockTitle(const char* title, bool checkFirstMatr) const
{
	if (checkFirstMatr) {
		if (!designParams()->firstMatr)
			return;

		designParams()->firstMatr = false;
	}

	char buffer[64];
	SPRINTF(buffer, " \n" BEG_OUT_BLOCK "%s: " END_OUT_BLOCK "\n", title);
	outString(buffer, outFile());
}

#if PRINT_SOLUTIONS
FClass2(CEnumerator, void)::printSolutions(const RowSolutionPntr pSolution, FILE* file, T nRow, bool markNextUsed, T nPartStart, T nPartEnd) const
{	
	if (!pSolution || nPartStart >= nPartEnd)
		return;

	const auto multiPortion = this->numParts() > 1;
	MUTEX_LOCK(out_mutex);
	for (auto i = nPartStart; i < nPartEnd; i++)
		(pSolution + i)->printSolutions(file, markNextUsed, nRow, i, multiPortion);

	MUTEX_UNLOCK(out_mutex);
}
#endif

#if USE_STRONG_CANONICITY_A
FClass2(CEnumerator, void)::checkUnusedSolutions(CRowSolution *pRowSolution)
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
	pRowSolution->resetSolutionIndex();
	CRowSolution *pRowSol = pRowSolution;
	size_t level = solIdx + 1;
	while (pRowSol->NextSolution(true) && pRowSol->solutionIndex() <= solIdx) {
		const auto *pCurrSol = pRowSol->currSolution();
		CColOrbit *pColOrb = MakeRow(pCurrSol, false);

		setCurrentRowNumb(nRow);
		setCurrUnforcedOrbPtr(nRow);
		setColOrbitCurr(pColOrb);
		canonChecker()->TestCanonicity(nRow, this, 0, &level, NULL, pRowSol);
		reset(nRow);
		setCurrentRowNumb(nRow - 1);
	}

	pRowSolution->setSolutionIndex(level - 1);
}
#endif

#if CANON_ON_GPU
FClass2(CEnumerator, size_t)::copyColOrbitInfo(S nRow) const
{
	// Function copy existing information about the orbits of columns  
	// into new structures, which could be used on GPU
	const auto pColOrbInfoBeg = CanonCheckerGPU()->ColOrbitData(t_CPU);
	auto pColOrbInfo = pColOrbInfoBeg;
	const auto colOrbit = colOrbits();
	const auto colOrbitIni = colOrbitsIni();
	while (nRow--) {
		const auto pColOrbitIni = colOrbitIni[nRow];
		auto pColOrbit = colOrbit[nRow];
		if (!pColOrbit) {
			*pColOrbInfo++ = ELEMENT_MAX;
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
