#pragma once
#include "SRGToolkit.h"
#include <cstring>

#pragma execution_character_set("utf-8")
class SrgSummary
{
public:
	SrgSummary(const std::string* logFolder, int iMode, const std::string* extFolder);
	void copyFiles(const std::string& inp, const std::string& out);
	void deleteFile(const std::string& tmp);
	void outSRG_info(int v, const SRGParam* graphParam, t_graphType graphType, int rank3, size_t grOrder, int srcGroupSize, int srcGroups, int srcAut);
private:
	std::string outFile;
	std::string outFile2;
	std::string tmpFile;
	std::string inpFile;
	int mode = 0;
};
