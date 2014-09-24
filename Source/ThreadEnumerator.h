
typedef enum {
	t_threadUndefined,
	t_threadLaunchFail,
	t_threadRunning,
	t_threadFinished,
	t_threadNotUsed
} t_threadCode;

class CEnumerator;
class CEnumInfo;

class CThreadEnumerator
{
public:
	CThreadEnumerator()								{ reset(); }
	~CThreadEnumerator()							{ release(); }
	void setupThreadForBIBD(const CEnumerator *pMaster, size_t nRow);
	void EnumerateBIBD(const CEnumerator *pMaster);
	inline t_threadCode code() const				{ return m_code; }
	inline CEnumerator *enumerator() const			{ return m_pEnum; }
	inline CEnumInfo *enumInfo() const				{ return m_pInfo; }
	inline void reInit()							{ release(); reset(); }
#if WRITE_MULTITHREAD_LOG
	inline void setThreadID(int id)					{ m_threadID = id; }
	inline int threadID() const					    { return m_threadID; }
#endif
	inline void setCode(t_threadCode code)			{ m_code = code; }
#if USE_BOOST && USE_POOL
	inline void setThread(boost::thread *pThread)	{ m_pTread = pThread; }
	inline boost::thread *getThread() const			{ return m_pTread; }
#endif
private:
	void release();
	void reset();

	CEnumInfo *m_pInfo;
	CEnumerator *m_pEnum;
	t_threadCode m_code;
#if USE_BOOST && USE_POOL
	boost::thread *m_pTread;
#endif
#if WRITE_MULTITHREAD_LOG
	int m_threadID;
#endif
};