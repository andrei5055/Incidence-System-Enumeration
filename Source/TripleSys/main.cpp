
#include <iostream>
#include "TripleSys.h"

int main()
{
	alldata sys(nPlayers);

#if nPlayers == 15
	const char* ivc = "";
#elif nPlayers == 21        
	const char* ivc = "";
#elif nPlayers == 27
	const char* ivc = ""
		/**/
 "     0   1   2    3   4   5    6   7   8    9  10  11   12  13  14   15  16  17   18  19  20   21  22  23   24  25  26 \n"
 "     0   3   6    1   4   7    2   5   8    9  12  15   10  13  16   11  14  19   17  21  24   18  22  25   20  23  26 \n"
 "     0   4   8    1   5   6    2   3   7    9  13  18   10  12  17   11  15  24   14  20  21   16  22  26   19  23  25 \n"
 "     0  10  26    1   9  22    2  14  18    3  15  19    4  12  24    5  20  25    6  17  23    7  11  16    8  13  21 \n"
		/**/
		;
#else
	const char* ivc = "";
#endif
	std::cout << "TripleSys3\n";

	sys.initStartValues(ivc);// can be used to start from previous result
	sys.Run();
	return 0;
}
/*
Result table:
 "     0   1   2    3   4   5    6   7   8    9  10  11   12  13  14   15  16  17   18  19  20 \n"
 "     0   5   8    1   9  14    2  10  18    3  11  19    4  12  15    6  16  20    7  13  17 \n"
 "     0   3   6    1  10  17    2   9  20    4  13  16    5  15  18    7  11  12    8  14  19 \n"
 "     0   4   7    1  15  19    2  11  17    3  10  16    5   9  13    6  14  18    8  12  20 \n"
 "     0   9  17    1   5   6    2  12  19    3  14  15    4  10  20    7  16  18    8  11  13 \n"
 "     0  10  14    1  12  16    2   5   7    3  13  20    4   9  19    6  11  15    8  17  18 \n"
 "     0  11  16    1   4   8    2  13  15    3   9  18    5  10  12    6  17  19    7  14  20 \n"
 "     0  12  18    1  11  20    2   3   8    4  14  17    5  16  19    6  10  13    7   9  15 \n"
 "     0  13  19    1   3   7    2  14  16    4  11  18    5  17  20    6   9  12    8  10  15 \n"
 "     0  15  20    1  13  18    2   4   6    3  12  17    5  11  14    7  10  19    8   9  16 \n"

Total time = 37360 ms
Total time = 32609 ms
TripleSys_5: Total time = 27938 ms
TripleSys_7: Total time = 30187 ms 
*/