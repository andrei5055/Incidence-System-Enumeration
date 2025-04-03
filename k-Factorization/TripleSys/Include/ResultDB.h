#pragma once

class ResultDB {
public:
	void UpdateEnumerationDB(char** pInfo, int len);
private:
	size_t getDirectory(char* dirName, size_t lenBuffer, bool rowNeeded) const;
	bool fileExists(const char* path, bool file) const;
	bool outFileIsValid(const struct stat& info, const char* pFileName = NULL) const { 
		return info.st_size > outFileValidSize(); 
	}
	size_t outFileValidSize() const					{ return 50; }
};
