
template<typename T>
void outMatrix(const T* c, int nl, int nc, int np, int ns, FILE* f, bool makeString, bool toScreen, int cntr=-1, const unsigned char* pDayPerm=NULL) {
	char buffer[512];
	const auto* endLine = makeString ? " \"\n" : "\n";
	for (int j = 0; j < nl; j++) {
		char* pBuf = buffer;
		SPRINTFD(pBuf, buffer, " \"");
		for (int i = 0; i < nc; i++, c++) {
			if (np && !(i % np))
				SPRINTFD(pBuf, buffer, " ");

			if (*c == -1 && !f)
				printfGreen(" %3d", *c);
			else
				SPRINTFD(pBuf, buffer, " %3d", *c);
		}

		if (cntr < 0) {
			if (j + 1 >= nl || ns <= 0 || ((j + 1) % ns) == 0)
				SPRINTFD(pBuf, buffer, endLine);
			else
				SPRINTFD(pBuf, buffer, " ");
		}
		else {
			if (cntr)
				SPRINTFD(pBuf, buffer, "\":  day =%2d\n", pDayPerm[j]);
			else
				SPRINTFD(pBuf, buffer, "\"\n");
		}

		_printf(f, toScreen, buffer);
	}
}

template<typename T>
class Table {
public:
	Table(char const *name, int nl, int nc, int ns = 0, int np = GroupSize, bool makeString = false, bool outCntr = false) :
		m_name(name), m_nl(nl), m_nc(nc), m_ns(ns), m_np(np), 
		m_makeString(makeString), m_bOutCntr(outCntr) {}
	void printTable(const T *c, bool outCntr = false, const char *fileName = NULL);
	inline void addCounterToTableName(bool val) { m_bOutCntr = val; }
private:
	const char *m_name;
	const int m_nl;
	const int m_nc;
	const int m_ns;
	const int m_np;
	const bool m_makeString;
	bool m_bOutCntr;          // When true, the counter will be added to the Table Name
public:
	int m_cntr = 0;
};

template<typename T>
void Table<T>::printTable(const T *c, bool outCntr, const char *fileName)
{
	char buffer[512], *pBuf = buffer;
	FOPEN_F(f, fileName, m_cntr ? "a" : "w");
	const auto* endLine = m_makeString ? " \"\n" : "\n";
	if (outCntr)
		m_cntr++;

	if (m_name && strlen(m_name) != 0) {
		if (outCntr && m_bOutCntr)
			SPRINTFD(pBuf, buffer, "%s %d:\n", m_name, m_cntr);
		else
			SPRINTFD(pBuf, buffer, "%s:\n", m_name);
	}

	_printf(f, true, buffer);
	outMatrix(c, m_nl, m_nc, m_np, m_ns, f, m_makeString, true);
	FCLOSE_F(f);
}

