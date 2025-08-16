#pragma once
#include "SRGToolkit.h"
#include <cstring>

#pragma execution_character_set("utf-8")
class SrgSummary
{
public:
	SrgSummary() {
		fopen_s(&f, "../../SRG.csv", "a");
		fseek(f, 0, SEEK_END);
		long initialSize = ftell(f);
		if (initialSize == 0)
			fprintf(f, "%s", "v,k,λ,μ,α,β,Src Groups,Src Group Size,Src Aut(M)\n");
	}
	~SrgSummary() {
		fclose(f);
	}
	void row(int v, int k, int λ, int μ, int α, int β, int srcGroups, int srcGroupSize, int srcAut) {
		if (v == -1)
			fprintf(f, ",");
		else
			fprintf(f, "%d,", v);
		fprintf(f, "%d,%d,%d,", k, λ, μ);
		if (α == -1)
			fprintf(f, "N/A,N/A,");
		else
			fprintf(f, "%d,%d,", α, β);
		fprintf(f, "%d,%d,%d\n", srcGroups, srcGroupSize, srcAut);
	}
private:
	FILE* f = NULL;
};
