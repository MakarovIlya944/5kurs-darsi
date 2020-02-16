#include "net.hpp"

void main()
{
	NetGenerator a = NetGenerator();
	a.isLoging = true;
	Net* net;
	Net v = a.GenerateFromFiles("D:\\Projects\\5kurs-darsi\\solver\\net.txt", "D:\\Projects\\5kurs-darsi\\solver\\border.txt", "D:\\Projects\\5kurs-darsi\\solver\\time.txt");
	int d = 4;
}