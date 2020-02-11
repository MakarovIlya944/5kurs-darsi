#include "MatrixGenerator.h"

double MatrixGenerator::CoefficientDiffusion(Pointd v)
{
	return 0;
	switch (curType)
	{
	case 1:
		return v.r / 2 + 3;
		break;
	case 2:
		return -v.r / 2 + 9;
		break;
	default:
		return 0;
		break;
	}
}

double MatrixGenerator::CoefficientSigma(Pointd v)
{
	return 1;
	switch (curType)
	{
	case 1:
		//return v.r + v.z;
		return 6 - v.r;
		break;
	case 2:
		return v.r - 6;
		break;
	default:
		return 0;
		break;
	}
}

double MatrixGenerator::RightPart(Pointd v,double t)
{
	return 6*t*t;
}

double MatrixGenerator::FirstCond(Pointd v)
{
	return u(v, GlobalTime[time]);
	switch (curType)
	{
	case 1:
		//return v.r*v.r + 1;
		//return v.r*v.r*v.r+2*v.r+4;
		return 1;
		break;
	case 2:
		//return v.z*v.z + 25;
		//return v.z*v.z*v.z*v.z+6*v.z+132;
		return v.z*v.z;
		break;
	case 3:
		//return v.r*v.r + 25;
		//return v.r*v.r*v.r+6*v.r+632;
		return 25;
		break;
	case 4:
		//return v.z*v.z + 1;
		//return v.z*v.z*v.z*v.z+2*v.z+4;
		return v.z*v.z;
		break;
	default:
		return 0;
		break;
	}
}

double MatrixGenerator::SecondCond(Pointd v)
{
	return -2 * GlobalTime[time]* GlobalTime[time];
	switch (curType)
	{
	case 1:
		return -2 * v.z;
		//return -2*(4*v.z*v.z*v.z+v.r+1);
		//return 5 * v.z*v.z*v.z*v.z;
		break;
	case 2:
		return 2 * v.r;
		//return 2 * (3*v.r * v.r + v.z + 1);
		//return 10 * v.r*v.r*v.r*v.r;
		break;
	case 3:
		//return 2*v.z;
		return 2 * (4 * v.z*v.z * v.z + v.r + 1);
		break;
	case 4:
		//return -2*v.r;
		return -2 * (3 * v.r * v.r + v.z + 1);
		break;
	default:
		return 0;
		break;
	}
}

double MatrixGenerator::ThirdCond(Pointd v)
{
	switch (curType)
	{
	case 1:
		//return -2*v.z;
		return v.r*v.r + 1 - v.z;
		break;
	case 2:
		//return 2*v.r;
		return 2 * (3 * v.r * v.r + v.z + 1);
		break;
	case 3:
		//return 2*v.z;
		return 2 * (4 * v.z*v.z * v.z + v.r + 1);
		break;
	case 4:
		//return -2*v.r;
		return -2 * (3 * v.r * v.r + v.z + 1);
		break;
	default:
		return 0;
		break;
	}
}

double MatrixGenerator::u(Pointd v, double t)
{
	return t*t*t;
}