#pragma once
#include "SRGToolkit.h"
#include <cstring>

#pragma execution_character_set("utf-8")
class SrgSummary
{
public:
	SrgSummary(std::string* logFolder, int iMode, std::string* extFolder);
	void copyFiles(std::string inp, std::string out);
	void deleteFile(std::string tmp);
	void outSRG_info(int v, const SRGParam* graphParam, t_graphType graphType, int rank3, size_t grOrder, int srcGroupSize, int srcGroups, int srcAut);
private:
	std::string outFile;
	std::string outFile2;
	std::string tmpFile;
	std::string inpFile;
	int mode = 0;
};
