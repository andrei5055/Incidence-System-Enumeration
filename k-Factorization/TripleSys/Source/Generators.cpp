#include "Table.h"

void reportNestedGroupCheckResult(int retVal, bool outToScreen) {
	if (retVal > 0) {
		printfRed("Nested groups check failed on row %d\n", retVal);
		myExit(1);
	}
	else {
		if (outToScreen && !retVal)
			printfGreen("Nested groups check OK\n");
	}
}

