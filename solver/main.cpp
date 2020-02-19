#ifdef _DEBUG
#include "net.hpp"

int main()
{
	NetGenerator a = NetGenerator();
	//a.isLoging = true;
	Net *v = a.GenerateFromFiles("D:\\Projects\\5kurs-darsi\\solver\\net.txt", "D:\\Projects\\5kurs-darsi\\solver\\border.txt", "D:\\Projects\\5kurs-darsi\\solver\\time.txt");
	Pointd tr = v->GlobalNet[1];
	int d = 4;
	return 0;
}
#endif 