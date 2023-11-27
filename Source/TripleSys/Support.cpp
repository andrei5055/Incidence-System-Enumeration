
#include <iostream>
#include "TripleSys.h"

#if UseBitmask
void checkbmask(alldata *s)
{
	for (int i = 0; i < nPlayers; i++)
	{
		if (s->selPlayers[i] != unset && (s->bmask & (1 << i)) != 0)
			i = i;
		else if (s->selPlayers[i] == unset && (s->bmask & (1 << i)) == 0)
			i = i;
	}
}
#endif
