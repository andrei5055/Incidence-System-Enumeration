#include <fstream>
#include <filesystem>
#include "TopGun.h"

void SizeParam::createFolderAndFileName(std::string& fn, const kSysParam* param, int tFolder, int nr, const std::string* fName) const
{
	namespace fs = std::filesystem;
	const char* folder = param->strVal[tFolder]->c_str();
	if (!folder || !strlen(folder))
		return;

	const char *uf, *fhdr = getFileNameAttr(param, &uf);
	fn = std::format("{}{}x{}x{}{}/", folder, m_numPlayers, nr, m_groupSize, uf);
	if (fs::create_directories(fn))
	{
		printfRed("*** Error: Cant create folder %s\n", fn.c_str());
		myExit(1);
	}
	if (fName)
		fn += fhdr + *fName;
}
