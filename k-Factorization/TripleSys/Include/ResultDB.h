#pragma once

class ResultDB {
public:
	void UpdateEnumerationDB(char** pInfo, int len);
private:
	size_t getDirectory(char* dirName, size_t lenBuffer, bool rowNeeded) const;
};
