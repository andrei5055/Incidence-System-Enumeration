#pragma once
#define nPlayers 9
#define GroupSize 3
#define nDays (nPlayers - 1) / (GroupSize - 1)
#define unset ((char)(-1))

struct alldata {
	int maxPlayers;
	int bmask;
	char maxResult[nDays][nPlayers];
	int maxDays;
	int maxDaysPlayers;
	int nLoops;
	bool noResult;
	char result[nDays][nPlayers];
	char links[nPlayers][nPlayers];
	char bklinks[nPlayers + 1][nPlayers][nPlayers];
	char devCounts[nPlayers];
	char selPlayers[nPlayers];
	char tmpPlayers[nPlayers];
	char indexPlayer[nPlayers];
	char preSet[nDays][nPlayers];
	int  iPlayer, iDay;

	void Init();
	bool Run(int numPlayers= nPlayers);
};

void adjust0(alldata* s);
bool setLinksAndDevCounts(alldata* s, char* p, int ip, char iset);
void printResult(alldata* s);
void init(alldata* s, const char** data);
void initPrevDayProcess(alldata* s);
void initCurrentDayProcess(alldata* s);
int check1(char result[][nPlayers], int ncolumns);
int compareMatrix(char result[][nPlayers], int ncolumns, char* transition);



