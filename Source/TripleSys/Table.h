
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
	int m_cntr = 0;
};

template<typename T>
void Table<T>::printTable(const T *c, bool outCntr, const char *fileName)
{
	char buffer[512], *pBuf = buffer;
	FILE* f = NULL;
	if (fileName && strlen(fileName))
		fopen_s(&f, fileName, m_cntr? "a" : "w");

	const auto* endLine = m_makeString ? " \\n\"\n" : "\n";
	if (outCntr)
		m_cntr++;

	if (m_name && strlen(m_name) != 0) {
		if (outCntr && m_bOutCntr)
			SPRINTF(pBuf, buffer, "%s %d:\n", m_name, m_cntr);
		else
			SPRINTF(pBuf, buffer, "%s:\n", m_name);
	}

	for (int j = 0; j < m_nl; j++)
	{
		if (m_makeString) SPRINTF(pBuf, buffer, " \" ");
		for (int i = 0; i < m_nc; i++, c++)
		{
			if (m_np > 0 && (i % m_np) == 0)
				SPRINTF(pBuf, buffer, " ");
			if (*c == -1 && !fileName)
				printfGreen(" %3d", *c);
			else
				SPRINTF(pBuf, buffer, " %3d", *c);
		}
		if (j + 1 >= m_nl || m_ns <= 0 || ((j + 1) % m_ns) == 0)
			SPRINTF(pBuf, buffer, endLine);
		else
			SPRINTF(pBuf, buffer, " ");

		_printf(f, true, pBuf = buffer);
	}
	FCLOSE(f);
}

