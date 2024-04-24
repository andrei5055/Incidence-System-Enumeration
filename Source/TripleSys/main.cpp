#include "TopGun.h"

int main()
{
	std::cout << "TripleSys 14.3\n";

	TopGun start(nPlayers, GroupSize);

	start.Run();

	std::cout << "\7" << std::endl; // play sound

	return 0;
}
