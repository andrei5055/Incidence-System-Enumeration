#include "TripleSys.h"

void convertLinksToResult(const char *ci, char* res, int np, int gs)
{
	const auto len = np * np;
	int nd = (np - 1) / (gs - 1);
	char* co = res;
	char *c = new char[len];
	memcpy(c, ci, len);
	char* a = new char[gs];
	for (int i = 0; i < nd; i++, co += np)
	{
		int m = 0;

		memset(co, 0, np);
		auto* cj = c;
		const auto gs_1 = gs - 1;
		for (int j = 0; j < np; j++, cj += np)
		{
			int cnt = 0;
			for (int k0 = 0; k0 < np; k0++)
			{
				if (cj[k0] == i)
				{
					if (cnt >= gs)
					{
						printf("too many players on day %d (%d, %d)\n", i, j, k0);
						break;
					}
					a[cnt++] = k0;
				}
			}
			if (cnt == gs_1)
			{
				if (m >= np - cnt)
				{
					printfRed("*** Error on day %d (%d)\n", i, m);
					break;
				}

				co[m] = j;
				memcpy(co + m + 1, a, gs_1);
				for (int i0 = 0; i0 < gs_1; i0++)
				{
					const auto idx0 = co[m + i0];
					if (idx0 == unset)
						continue;

					for (int i1 = i0 + 1; i1 < gs; i1++)
					{
						const auto idx1 = co[m + i1];
						if (idx1 != unset) {
							c[idx0 * np + idx1] = unset;
							c[idx1 * np + idx0] = unset;
						}
					}
				}
				m += gs;
			}
		}
	}
	delete[] a;
	delete[] c;
}