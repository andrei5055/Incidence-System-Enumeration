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
		const auto initialSize = ftell(f);
		if (!initialSize)
			fprintf(f, "%s", "v,k,λ,μ,α,β,Rank 3,k_c,λ_c,μ_c,Aut(G),Group Size,Number of Groups,Src Aut(M)\n");
	}
	~SrgSummary() {
		fclose(f);
	}
	void outSRG_info(int v, const SRGParam* graphParam, t_graphType graphType, int rank3, size_t grOrder, int srcGroupSize, int srcGroups, int srcAut) {
		char buf[1024], * pBuf = buf;
		SPRINTFD(pBuf, buf, "%d,", v);
//		fprintf(f, "%d,", v);

		const auto k = graphParam->k;
		const auto λ = graphParam->λ;
		const auto μ = graphParam->μ;
		SPRINTFD(pBuf, buf, "%d,%d,%d,", k, λ, μ);
//		fprintf(f, "%d,%d,%d,", k, λ, μ);
		if (graphType == t_4_vert)
			SPRINTFD(pBuf, buf, "%d,%d", graphParam->α, graphParam->β);
		//			fprintf(f, "%d,%d,", graphParam->α, graphParam->β);
		else
			SPRINTFD(pBuf, buf, "n/a,n/a");
//			fprintf(f, "n/a,n/a,");

		const auto rank2Symb = rank3 == 1 ? '+' : (rank3? '-' : ' ');
		const auto v_2k = v - 2 * k;
		const auto k_c = v - k - 1;
		const auto λ_c = v_2k + μ - 2;
		const auto μ_c = v_2k + λ;
		SPRINTFD(pBuf, buf, ",%c,%d,%d,%d,", rank2Symb, k_c, λ_c, μ_c);
		if (grOrder != -1)
			SPRINTFD(pBuf, buf, "%zu", grOrder);

		SPRINTFD(pBuf, buf, ",%d,%d,%d\n", srcGroupSize, srcGroups, srcAut);
		fprintf(f, buf);
//		fprintf(f, "%c,%d,%d,%d,%d,%d,%d\n", rank2Symb, k_c, λ_c, μ_c, srcGroupSize, srcGroups, srcAut);
	}
private:
	FILE* f = NULL;
};
