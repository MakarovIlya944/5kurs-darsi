#include "Interface.h"

void Interface::InitFirstTime()
{
	for (int i(0); i < nx; i++)
		Q[0][i] = u(NetGen.GlobalNet[i],NetGen.GlobalTime[0]);

	for (int i(0); i < nx; i++)
		Q[1][i] = u(NetGen.GlobalNet[i], NetGen.GlobalTime[1]);
}

double Interface::u(Pointd v, double t)
{
	return t*t*t;
}