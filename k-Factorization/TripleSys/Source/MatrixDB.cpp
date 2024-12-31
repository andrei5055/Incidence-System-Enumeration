#include "MatrixDB.h"
#include <ctime>

MatrDescr::MatrDescr(size_t groupOrder, const char *cyclesDescr, size_t numMatrix) : m_groupOrder(groupOrder) {
	m_numMatrix = numMatrix;
	const auto len = strlen(cyclesDescr) + 1;
	m_cyclesDescr = new char[len];
	memcpy(m_cyclesDescr, cyclesDescr, len);
}

MatrixDB::~MatrixDB() {
	if (descrStorage()) {
		auto& storage = *descrStorage();
		for (auto it = begin(storage); it != end(storage); ++it)
			delete *it;

		delete descrStorage();
	}
}

void MatrixDB::transferMatrixDB(MatrixDB* matrDB) {
	delete descrStorage();
	setDescrStorage(matrDB->descrStorage());
	matrDB->setDescrStorage(NULL);
}

void MatrixDB::addMatrixDB(const MatrixDB* matrDB) {
	//auto& storage = *descrStorage();
	//auto it = begin(storage);
	if (matrDB->descrStorage()) {
		for (const auto& record : *matrDB->descrStorage()) {
			addMatrix(record->groupOrder(), record->cycleDescr(), record->numMatrix());
		}
	}
}

void MatrixDB::addMatrix(size_t groupOrder, const char* cyclesDescr, size_t numMatrix) {
	bool addNew = true;
	if (descrStorage()) {
		auto& storage = *descrStorage();
		for (auto it = begin(storage); it != end(storage); ++it) {
			if ((*it)->groupOrder() < groupOrder)
				continue;

			int cmp = 0;
			if ((*it)->groupOrder() > groupOrder || (cmp = strcmp((*it)->cycleDescr(), cyclesDescr)) != 0) {
				if (cmp < 0)
					continue;

				storage.insert(it, new MatrDescr(groupOrder, cyclesDescr, numMatrix));
			}
			else
				(*it)->addNumMatrix(numMatrix);

			addNew = false;
			break;
		}
	}
	else
		setDescrStorage(new std::vector<MatrDescr*>);

	if (addNew)
		descrStorage()->push_back(new MatrDescr(groupOrder, cyclesDescr, numMatrix));
}

void out_date_time(FILE* f) {
	time_t timestamp = time(NULL);
	struct tm datetime;
	localtime_s(&datetime, &timestamp);

	char output[50];
	strftime(output, sizeof(output), "%m/%d/%y: %H:%M:%S", &datetime);
	fprintf(f, "\n%s %s\n", DATE_TIME_TAG, output);
}

void MatrixDB::reportResult(FILE* f) const {
	out_date_time(f);
	if (descrStorage()) {
		char line[61];
		memset(line, '-', sizeof(line));
		line[sizeof(line) - 1] = '\0';
#define My_FPRINTF(f, format, ...) fprintf(f, "     "##format, __VA_ARGS__)
		My_FPRINTF(f, "%s\n", line);
		My_FPRINTF(f, "  |Aut(M)|:              Cycles:            # of matrices:\n");
		My_FPRINTF(f, "%s\n", line);
		size_t total = 0;
		for (const auto& record : *descrStorage()) {
			My_FPRINTF(f, " %6zd  %34s  %8zd\n", record->groupOrder(), record->cycleDescr(), record->numMatrix());
			total += record->numMatrix();
		}

		My_FPRINTF(f, "%s\n", line);
		My_FPRINTF(f, "%42s   %8zd\n", "Total:", total);
		My_FPRINTF(f, "%s\n\n\n", line);
	}
	else {
		My_FPRINTF(f, " *** Resulting matrix DB is empty *** \n");
	}
}