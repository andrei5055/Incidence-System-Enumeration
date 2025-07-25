#pragma once
#include <string> 
#include <vector>
#include <mutex>
#include <memory>
#include "Global.h"

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
	inline void setDescrStorage(std::vector<T*>* ptr) {
		m_descrStorage.reset(ptr);
	}
	std::vector<std::unique_ptr<T>>& descrStorage() {
		return m_descrStorage;
	}

	const std::vector<std::unique_ptr<T>>& descrStorage() const {
		return m_descrStorage;
	}

	void addObjDescriptor(size_t groupOrder, const char *objDescr, size_t numMatrix = 1, size_t sourceID = 0);
	size_t reportResult(FILE* f, bool writeTimeStamp = true) const;
	inline void setTableTitle(const char* title) { m_tableTile = title; }
protected:
	virtual void tableTitle(FILE* f) const { My_FPRINTF(f, "%s\n", m_tableTile.c_str()); }
	virtual const char* columnTitles() const = 0;
	virtual int tableWidth() const = 0;
	virtual size_t printRecord(FILE *f, const ObjDescr *record) const = 0;
	virtual bool outTotal() const = 0;
	virtual void printTotal(FILE* f, size_t total) const = 0;
	virtual void reportEmptyDB(FILE* f) const = 0;
	inline auto totalObjects() const		{ return m_totalObjects; }
private:
	//std::unique_ptr<std::vector<T*>> m_descrStorage;
	std::vector<std::unique_ptr<T>> m_descrStorage;
	size_t m_totalObjects = 0;
	std::string m_tableTile;
};

template <class T>
void DescrStorage<T>::addObjDescriptor(size_t groupOrder, const char* objDescr, size_t numMatrix, size_t sourceID) {
	this->m_totalObjects++;
	for (size_t i = 0; i < m_descrStorage.size(); ++i) {
		T* record = m_descrStorage[i].get();

		const auto currGroupOrder = record->groupOrder();
		if (currGroupOrder < groupOrder)
			continue;

		int cmp = 0;
		if (currGroupOrder > groupOrder || (cmp = strcmp(record->objDescr(), objDescr)) != 0) {
			if (cmp < 0)
				continue;

			m_descrStorage.insert(m_descrStorage.begin() + i,
				std::make_unique<T>(groupOrder, objDescr, numMatrix, sourceID));
			return;
		}
		else {
			record->incNumObj(numMatrix);
			record->addSourceObj(sourceID);
			return;
		}
	}

	// Insert at end if not found
	m_descrStorage.emplace_back(std::make_unique<T>(groupOrder, objDescr, numMatrix, sourceID));
}

template <class T>
size_t DescrStorage<T>::reportResult(FILE* f, bool writeTimeStamp) const {
	void out_date_time(FILE * f);
	if (writeTimeStamp)
		out_date_time(f);

	size_t total = 0;
	const auto& storage = descrStorage();
	if (!storage.empty()) {
		tableTitle(f);
		const auto lenLine = tableWidth();
		auto *line = new char[lenLine];
		memset(line, '-', lenLine - 1);
		line[lenLine - 1] = '\0';
		My_FPRINTF(f, "%s\n", line);
		My_FPRINTF(f, "%s\n", columnTitles());
		My_FPRINTF(f, "%s\n", line);

		for (const auto& record : storage) {
			record->addSourceObj(-1);
			total += printRecord(f, record.get());
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
	void updateMatrixDB(MatrixDB* db);
	void transferMatrixDB(MatrixDB* matrDB);
	void addMatrixDB(MatrixDB* matrDB);
protected:
	virtual int tableWidth() const			{ return 61; }	
	virtual const char* columnTitles() const{ return "  |Aut(M)|:  # of matrices:          Cycles:"; }
	virtual size_t printRecord(FILE* f, const ObjDescr* record) const {
		My_FPRINTF(f, "%8zd      %8zd        %-34s\n", record->groupOrder(), record->numObjects(), record->objDescr());
		return record->numObjects();
	}
	virtual bool outTotal() const			{ return true; }
	virtual void printTotal(FILE *f, size_t total) const {
		My_FPRINTF(f, "%11s   %8zd\n", "Total:", total);
	}
	virtual void reportEmptyDB(FILE *f) const {
		My_FPRINTF(f, " *** Resulting matrix DB is empty *** \n");
	}
private:
	mutable std::mutex m_descrMutex[2];
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
	inline void setGraphType(int graphType)		{ m_graphType = graphType; }
	inline void setGraphType(t_graphType type)	{ m_grType = type; }
protected:
	virtual int tableWidth() const				{ return 65; }
	virtual const char* columnTitles() const	{
		return "  |Aut(G)|:    # of graphs:    Graph type:        Sources:"; }
	virtual size_t printRecord(FILE* f, const ObjDescr* record) const {
		My_FPRINTF(f, "%10zd         %5zd      %-18s   %-23s\n", record->groupOrder(), record->numObjects(), record->objDescr(), record->getSources());
		return record->numObjects();
	}
	virtual bool outTotal() const				{ return true; }
	virtual void printTotal(FILE* f, size_t total) const {
		My_FPRINTF(f, "%18s %5zd                      %9zd\n", "Total:", total, totalObjects());
	}
	virtual void reportEmptyDB(FILE* f) const	{
		My_FPRINTF(f, " *** The graphs of type %d are %s and therefore have not been constructed. *** \n\n\n", 
			m_graphType, (m_grType == t_regular? "merely regular" : "complete"));
	}
private:
	int m_graphType = 0;
	t_graphType m_grType;
};
