#include <string.h>
#include <stdio.h>

#include "ResultDB.h"

#define countof(x) sizeof(x)/sizeof(x[0])

void ResultDB::UpdateEnumerationDB(char** pInfo, int len)
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

	FILE* dbFile = NULL;
	int i = -1;
	while (!dbFile && ++i < 2) {
		FOPEN(file, enumerationDB, i ? "w" : "r");
		dbFile = file;
	}

	if (!dbFile)
		return;

	char key[32], adjustedKey[64];// , keyCmp[64];
	this->getEnumerationObjectKey(key, countof(key));
	// the lexicographical order of key's could be adjusted for some type of designs using:
	const char* pAdjKey = this->getEnumerationObjectKeyA(adjustedKey, countof(adjustedKey) / 2);

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
	const auto lenKeyAdj = pAdjKey ? strlen(pAdjKey) : 0;
	bool firstLine = true;   // We don't need to compare key with the first line
	while (fgets(buffer, countof(buffer), dbFile)) {
		if (buffer[0] != ';') {
			// Not a comment line
			if (compareFlag && !firstLine) {
				if (buffer[0] && !strchr("_-", buffer[0]))
					resCmp = pAdjKey ? compareEnumerationDB_record(buffer) : strncmp(buffer, key, lenKey);
				else
					resCmp = 1;   // EOF found

				if (resCmp >= 0) {
					const char* pComment = !resCmp ? strstr(buffer, " >> ") : NULL;
					outKeyInfo(key, pInfo, f, pComment);
					compareFlag = false;
					if (!resCmp)
						continue;
				}
			}
		}

		firstLine = false;
		if (f && fputs(buffer, f) < 0) {
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

size_t ResultDB::getDirectory(char* dirName, size_t lenBuffer, bool rowNeeded) const
{
	const auto pParam = designParams();
	lenBuffer--;		// Reserving 1 byte for last '/'

	auto len = SNPRINTF(dirName, lenBuffer, "%s", pParam->strParam[t_workingDir].c_str());
	SET_DIRECTORY(dirName);

	const auto* pDirName = this->getTopLevelDirName();
	if (pDirName) {
		len += SNPRINTF(dirName + len, lenBuffer - len, "%s/", pDirName);
		SET_DIRECTORY(dirName);
	}

	if (rowNeeded) {
		auto rowNumb = getInSys()->rowNumbExt();
		if (pParam->objType == t_objectType::t_SemiSymmetricGraph)
			rowNumb *= pParam->r / pParam->k;

		len += SNPRINTF(dirName + len, lenBuffer - len, "V =%4d/", rowNumb);
		SET_DIRECTORY(dirName);
	}

	return len;
}