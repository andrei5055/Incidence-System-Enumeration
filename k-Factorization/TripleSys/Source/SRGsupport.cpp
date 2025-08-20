#pragma once
#include <fstream>
#include <iostream> 

#include <filesystem>
#include <cstring>
#include <SRGsupport.h>
using namespace std;

#pragma execution_character_set("utf-8")

SrgSummary::SrgSummary(string* logFolder, int iMode, string* extFile) {
	mode = iMode;
	if (!mode)
		return;
	inpFile = mode == 1 ? (extFile ? *extFile : *logFolder + "/SRG.csv") : "";
	tmpFile = *logFolder + "/SRGtmp.csv";
	outFile = *logFolder + "/SRG.csv";
	outFile2 = extFile ? *extFile : *logFolder + "/SRG.csv";
}

void SrgSummary::outSRG_info(int v, const SRGParam* graphParam, t_graphType graphType, int rank3, size_t grOrder, int srcGroupSize, int srcGroups, int srcAut) {
	if (!mode)
		return;

	FILE* fInp = 0, * fTmp = 0, * fOut = 0, * fOut2 = 0;
	if (mode == 1) {
		fopen_s(&fInp, inpFile.c_str(), "r");
		if (!fInp)
			inpFile = "";
	}
	fopen_s(&fTmp, tmpFile.c_str(), "w+");
	if (!fTmp) {
		printf("*** Error: Can't open temporary file to save SRG summary (%s)\n", tmpFile.c_str());
		exit(1);
	}

	char buf[1024], bufInp[1024];
	if (fInp)
		fgets(bufInp, sizeof(bufInp), fInp);

	fprintf(fTmp, "%s", "v,k,λ,μ,α,β,Rank 3,kc,λc,μc,Aut(G),Group Size,Number of Groups,Src Aut(M)\n");

	const auto k = graphParam->k;
	const auto λ = graphParam->λ;
	const auto μ = graphParam->μ;
	// form new row in 'buf'
	char* pBuf = buf;
	SPRINTFD(pBuf, buf, "%d,%d,%d,%d,", v, k, λ, μ);

	if (graphType == t_4_vert)
		SPRINTFD(pBuf, buf, "%d,%d", graphParam->α, graphParam->β);
	else
		SPRINTFD(pBuf, buf, "n/a,n/a");

	const auto rank2Symb = rank3 == 1 ? '+' : (rank3 ? '-' : ' ');
	const auto v_2k = v - 2 * k;
	SPRINTFD(pBuf, buf, ",%c,%d,%d,%d,", rank2Symb, v - k - 1, v_2k + μ - 2, v_2k + λ);
	if (grOrder != -1)
		SPRINTFD(pBuf, buf, "%zu", grOrder);
	SPRINTFD(pBuf, buf, ",%d,%d,%d\n", srcGroupSize, srcGroups, srcAut);

	// merge (if requested) new row with existing SRG.csv
	bool bNewRowRecorded = false;
	int vOld = -1, kOld = -1;
	while (1) {
		char* pStr = fInp ? fgets(bufInp, sizeof(bufInp), fInp) : NULL;
		if (!bNewRowRecorded) {
			if (pStr && strcmp(buf, bufInp) == NULL) {
				// new row is the same as current row from SRG.csv, ignore new row
				bNewRowRecorded = true;
			}
			else {
				if (pStr)
					sscanf_s(bufInp, "%d,%d", &vOld, &kOld);
				if (!pStr || vOld > v || (vOld == v && kOld > k)) {
					fprintf(fTmp, buf);
					bNewRowRecorded = true;
				}
			}
		}
		if (!pStr)
			break;
		fputs(bufInp, fTmp);
	}

	FCLOSE_F(fInp);
	FCLOSE_F(fTmp);

	copyFiles(tmpFile, outFile);
	copyFiles(tmpFile, outFile2);
	deleteFile(tmpFile);
}
void SrgSummary::deleteFile(string tmp)
{
	if (tmp.length()) {
		try {
			if (filesystem::remove(tmp)) {
				cout << "File deleted successfully: " << tmp << endl;
			}
			else {
				cout << "File not found or could not be deleted: " << tmp << endl;
			}
		}
		catch (const filesystem::filesystem_error& e) {
			cerr << "Error deleting file: " << e.what() << endl;
		}
	}
}
void SrgSummary::copyFiles(string inp, string out)
{
	if (inp.length() && out.length()) {
		filesystem::path source_path = inp;
		filesystem::path destination_path = out;

		try {
			filesystem::copy(source_path, destination_path, filesystem::copy_options::overwrite_existing);
			cout << "File copied successfully from " << source_path << " to " << destination_path << endl;
		}
		catch (const filesystem::filesystem_error& e) {
			cerr << "Error copying file: " << e.what() << endl;
		}
	}
}
