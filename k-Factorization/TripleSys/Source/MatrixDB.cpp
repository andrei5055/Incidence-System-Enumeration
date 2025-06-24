#include "MatrixDB.h"
#include <ctime>

void MatrixDB::transferMatrixDB(MatrixDB* matrDB) {
	delete descrStorage();
	setDescrStorage(matrDB->descrStorage());
	matrDB->setDescrStorage(NULL);
}

void MatrixDB::addMatrixDB(const MatrixDB* matrDB) {
	if (matrDB->descrStorage()) {
		for (const auto& record : *matrDB->descrStorage()) {
			addObjDescriptor(record->groupOrder(), record->objDescr(), record->numObjects());
		}
	}
}

void out_date_time(FILE* f) {
	time_t timestamp = time(NULL);
	struct tm datetime;
	localtime_s(&datetime, &timestamp);

	char output[50];
	strftime(output, sizeof(output), "%m/%d/%y: %H:%M:%S", &datetime);
	fprintf(f, "\n%s %s\n", DATE_TIME_TAG, output);
}

void GraphDescr::addSourceObj(size_t sourceID) {
	if (sourceID != m_prevID + 1) {
		if (m_prevID == -1)
			m_source = std::to_string(sourceID);
		else {
			if (m_cntr > 1)
				m_source += "-" + std::to_string(m_prevID);
			else
				if (m_cntr == 1 && sourceID == -1)
					sourceID = m_prevID;

			if (sourceID != -1)
				m_source += "," + std::to_string(sourceID);
		}

		m_cntr = 0;
	}
	else
		m_cntr++;

	m_prevID = sourceID;
}