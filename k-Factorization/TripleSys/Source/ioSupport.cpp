#include <fstream>
#include <filesystem>
#include "TopGun.h"

void SizeParam::createFolderAndFileName(std::string& fn, const kSysParam* param, int tFolder, int nr, const std::string& fName) const
{
	namespace fs = std::filesystem;
	const char* folder = param->strVal[tFolder]->c_str();
	if (!folder || !strlen(folder))
		return;

	const char* uf, * fhdr = getFileNameAttr(param, &uf);
	char buffer[32];
	if (param->completeGraph())
		snprintf(buffer, sizeof(buffer), "Complete");
	else
		snprintf(buffer, sizeof(buffer), "%d-Partite", param->partiteNumb());

	fn = std::format("{}{}_graphs/{}/{}x{}x{}{}/", folder, buffer, param->partitionSize(), m_numPlayers, nr, m_groupSize, uf);
	if (!fs::exists(fn)) {
		try {
			if (!fs::create_directories(fn)) {
				const auto err = GetLastError();
				if (err != ERROR_ALREADY_EXISTS) {
					printfRed("*** Error: Unable to create folder \"%s\" (system error %d(0x%x)\n", fn.c_str(), err, err);
					myExit(err ? err : 1);
				}
			}
		}
		catch (...) {
			printfRed("*** Error: Unable to create folder \"%s\"\n", fn.c_str());
			myExit(1);
		}
	}
	if (!fName.empty())
		fn += fhdr + fName;
}
