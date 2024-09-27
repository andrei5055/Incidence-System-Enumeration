#include <iostream>
#include "TripleSys.h"
#include <windows.h>
#include <conio.h>
void setConsoleOutputMode()
{
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	HANDLE hConsoleIn = GetStdHandle(STD_INPUT_HANDLE);
	INPUT_RECORD irInBuf[128];
	DWORD cNumRead = 0, fdOldMode = 0, i;
	static bool bFixedOut = false;
	//if (_kbhit())
	//if (GetConsoleMode(hConsoleIn, &fdOldMode))
	{
		SetConsoleMode(hConsoleIn, fdOldMode | ENABLE_WINDOW_INPUT);
		GetNumberOfConsoleInputEvents(hConsoleIn, &cNumRead);

		if (cNumRead > 0 && ReadConsoleInput(hConsoleIn, irInBuf, 128, &cNumRead) && cNumRead > 0)
		{
			for (i = 0; i < cNumRead; i++)
			{
				if (irInBuf[i].EventType == KEY_EVENT)
				{
					bFixedOut = !bFixedOut;
					break;
				}
			}
		}
	}
	if (bFixedOut)
	{
		COORD pos = { 0, 6 };
		SetConsoleCursorPosition(hConsole, pos);
		printfGreen("Press any button to stop/start fixed output\n");
	}
}