#pragma once
#include <string> 
#include <vector>

#define DATE_TIME_TAG   "===>>"
#define My_FPRINTF(f, format, ...) fprintf(f, "     "##format, __VA_ARGS__)

class ObjDescr {
public:
	ObjDescr(size_t groupOrder, const char* objDescr, size_t numObj, size_t sourceID) : m_groupOrder(groupOrder) {
		m_numObj = numObj;
		const auto len = strlen(objDescr) + 1;
		m_objDescr = new char[len];
		memcpy(m_objDescr, objDescr, len);
	}
	~ObjDescr()								{ delete[] objDescr(); }
	inline auto groupOrder() const			{ return m_groupOrder; }
	inline auto numObjects() const			{ return m_numObj; }
	inline void incNumObj(size_t n)			{ m_numObj += n; }
	inline char* objDescr() const			{ return m_objDescr; }
	virtual const char* getSources() const	{ return NULL; }
	virtual void addSourceObj(size_t sourceID) {}
private:
	const size_t m_groupOrder;	// order of automorphism group
	size_t m_numObj;			// number of objects with a given groupOrder/objDescr
	char* m_objDescr = NULL;	// obj description
};

class MatrDescr : public ObjDescr {
public:
	MatrDescr(size_t groupOrder, const char *cyclesDescr, size_t numMatrix, size_t sourceID) :
		ObjDescr(groupOrder, cyclesDescr, numMatrix, sourceID) {}
};

template <class T>
class DescrStorage {
public:
	~DescrStorage();
	inline void setDescrStorage(std::vector<T*>* pntr) {
		m_descrStorage = pntr;
	}
	inline auto descrStorage() const { return m_descrStorage; }
	void addObjDescriptor(size_t groupOrder, const char *objDescr, size_t numMatrix = 1, size_t sourceID = 0);
	size_t reportResult(FILE* f, bool writeTimeStamp = true) const;
protected:
	virtual const char* tableTitle() const = 0;
	virtual int tableWidth() const = 0;
	virtual size_t printRecord(FILE *f, const ObjDescr *record) const = 0;
	virtual bool outTotal() const = 0;
	virtual void printTotal(FILE* f, size_t total) const = 0;
	virtual void reportEmptyDB(FILE* f) const = 0;
private:
	std::vector<T*>* m_descrStorage = NULL;
};

template <class T>
DescrStorage<T>::~DescrStorage() {
	if (descrStorage()) {
		auto& storage = *descrStorage();
		for (auto it = begin(storage); it != end(storage); ++it)
			delete* it;

		delete descrStorage();
	}
}

template <class T>
void DescrStorage<T>::addObjDescriptor(size_t groupOrder, const char* objDescr, size_t numMatrix, size_t sourceID) {
	if (descrStorage()) {
		auto& storage = *descrStorage();
		for (auto it = begin(storage); it != end(storage); ++it) {
			if ((*it)->groupOrder() < groupOrder)
				continue;

			int cmp = 0;
			if ((*it)->groupOrder() > groupOrder || (cmp = strcmp((*it)->objDescr(), objDescr)) != 0) {
				if (cmp < 0)
					continue;

				storage.insert(it, new T(groupOrder, objDescr, numMatrix, sourceID));
			}
			else {
				(*it)->incNumObj(numMatrix);
				(*it)->addSourceObj(sourceID);
			}
			return;
		}
	}
	else
		setDescrStorage(new std::vector<T *>);

	descrStorage()->push_back(new T(groupOrder, objDescr, numMatrix, sourceID));
}

template <class T>
size_t DescrStorage<T>::reportResult(FILE* f, bool writeTimeStamp) const {
	void out_date_time(FILE * f);
	if (writeTimeStamp)
		out_date_time(f);

	size_t total = 0;
	if (descrStorage()) {
		const auto lenLine = tableWidth();
		auto *line = new char[lenLine];
		memset(line, '-', lenLine - 1);
		line[lenLine - 1] = '\0';
		My_FPRINTF(f, "%s\n", line);
		My_FPRINTF(f, "%s\n", tableTitle());
		My_FPRINTF(f, "%s\n", line);

		for (const auto& record : *descrStorage()) {
			record->addSourceObj(-1);
			total += printRecord(f, record);
		}

		if (outTotal()) {
			My_FPRINTF(f, "%s\n", line);
			printTotal(f, total);
		}
		My_FPRINTF(f, "%s\n\n\n", line);
		delete[] line;
	}
	else {
		reportEmptyDB(f);
	}
	return total;
}
class MatrixDB : public DescrStorage <MatrDescr> {
public:
	void transferMatrixDB(MatrixDB* matrDB);
	void addMatrixDB(const MatrixDB* matrDB);
protected:
	virtual int tableWidth() const			{ return 61; }	
	virtual const char* tableTitle() const	{ return "  |Aut(M)|:              Cycles:            # of matrices:"; }
	virtual size_t printRecord(FILE* f, const ObjDescr* record) const {
		My_FPRINTF(f, " %6zd  %34s  %8zd\n", record->groupOrder(), record->objDescr(), record->numObjects());
		return record->numObjects();
	}
	virtual bool outTotal() const			{ return true; }
	virtual void printTotal(FILE *f, size_t total) const {
		My_FPRINTF(f, "%42s   %8zd\n", "Total:", total);
	}
	virtual void reportEmptyDB(FILE *f) const {
		My_FPRINTF(f, " *** Resulting matrix DB is empty *** \n");
	}
};

class GraphDescr : public ObjDescr {
public:
	GraphDescr(size_t groupOrder, const char* graphDescr, size_t numObj, size_t sourceID) : ObjDescr(groupOrder, graphDescr, numObj, sourceID) {
		addSourceObj(sourceID);
	}
	virtual void addSourceObj(size_t sourceID);
protected:
	inline const char *getSources() const	{ return m_source.c_str(); }
private:
	std::string m_source;
	size_t m_prevID = (size_t)(-1);
	int m_cntr = 0;
};

class GraphDB : public  DescrStorage <GraphDescr> {
public:
	inline void setGraphType(int graphType) { m_graphType = graphType; }
protected:
	virtual int tableWidth() const			{ return 61; }
	virtual const char* tableTitle() const	{
		return "  |Aut(G)|:    # of graphs:    Graph type:       Sources:"; }
	virtual size_t printRecord(FILE* f, const ObjDescr* record) const {
		My_FPRINTF(f, " %6zd            %5zd     %-22s %-23s\n", record->groupOrder(), record->numObjects(), record->objDescr(), record->getSources());
		return record->numObjects();
	}
	virtual bool outTotal() const			{ return true; }
	virtual void printTotal(FILE* f, size_t total) const {
		My_FPRINTF(f, "%18s %5zd\n", "Total:", total);
	}
	virtual void reportEmptyDB(FILE* f) const {
		My_FPRINTF(f, " *** The graphs of type %d are merely regular and therefore have not been constructed. *** \n\n\n", m_graphType);
	}
private:
	int m_graphType = 0;
};
