#pragma once
#include <Windows.h>
#ifdef _DEBUG
#else
#include <pybind11/pybind11.h>
#endif 
#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <iterator>
#include <tuple>

using namespace std;

struct Pointd
{
	double r;
	double z;

	Pointd()
	{
		r = 0;
		z = 0;
	}

	Pointd(double R, double Z)
	{
		r = R;
		z = Z;
	}

	Pointd operator +(Pointd a)
	{
		return Pointd(r + a.r, z + a.z);
	}

	Pointd operator -(Pointd a)
	{
		return Pointd(r - a.r, z - a.z);
	}

};

struct Pointi
{
	int r;
	int z;
};

struct Net {
public:
	//rNum - num of points at r; zNum - num of points at z
	int rNum, zNum;
	//Num - num of points
	int Num;
	//NumT - num of time layers
	int NumT;
	//NumEl - num of elements
	int NumEl;
	//GlobalNet - koords of all points (r,z)
	Pointd* GlobalNet;
	//time net
	double* GlobalTime;
	//nvtr - elements
	int* nvtr;
	//nvcat - materials
	int* nvcat;
	//nvr1 - first border condition
	list<int> nvr1;
	//nvr2 - second border condition
	list<int> nvr2;
	//nvr3 - third border condition
	list<int> nvr3;

	~Net();
};

class NetGenerator
{
public:

	bool isLoging = false;

	//isSavingNet - if true save data fo file
	bool isSavingNet = false;
	const char* fileName = "default.txt";

	Net Generate(vector<double> R, vector<double> Z, double step, vector<double> r_steps, vector<double> z_steps, vector<int> materialIds, vector<Pointi> materials, vector<int> B, vector<int> V, vector<Pointi> border, vector<double> T, double t_step, vector<double> t_steps);

	//FNameN - net.txt , FNameB - border.txt , FNameT - time.txt
	Net* GenerateFromFiles(const char* FNameN, const char* FNameB, const char* FNameT);

	void GenerationNVTR(Net &net);
};

#ifdef _DEBUG

#else

namespace py = pybind11;

PYBIND11_MODULE(cpp_solver, m) {
	py::class_<NetGenerator>(m, "NetGenerator")
		.def(py::init<>())
		.def("GenerateFromFiles", &NetGenerator::GenerateFromFiles, R"pbdoc(FNameN - net.txt, FNameB - border.txt, FNameT - time.txt)pbdoc")
		.def("Generate", &NetGenerator::Generate, R"pbdoc(
Points at r
Points at z
Default step
Step and koefficient per elemens at r
Step and koefficient per elemens at z
Materials border
)pbdoc")
.def_readwrite("isLoging", &NetGenerator::isLoging);

	py::class_<Net>(m, "Net")
		.def(py::init<>())
		.def_readwrite("rNum", &Net::rNum, R"pbdoc(rNum - num of points at r)pbdoc")
		.def_readwrite("zNum", &Net::zNum, R"pbdoc(zNum - num of points at z)pbdoc")
		.def_readwrite("Num", &Net::Num, R"pbdoc(Num - num of points)pbdoc")
		.def_readwrite("NumT", &Net::NumT, R"pbdoc(NumT - num of time layers)pbdoc")
		.def_readwrite("NumEl", &Net::NumEl, R"pbdoc(NumEl - num of elements)pbdoc")
		.def_readwrite("GlobalNet", &Net::GlobalNet, R"pbdoc(GlobalNet - koords of all points (r,z))pbdoc")
		.def_readwrite("GlobalTime", &Net::GlobalTime, R"pbdoc(GlobalTime - all time layers)pbdoc")
		.def_readwrite("nvtr", &Net::nvtr, R"pbdoc(nvtr - all elements)pbdoc")
		.def_readwrite("nvcat", &Net::nvcat, R"pbdoc(nvcat - all materials)pbdoc")
		.def_readwrite("nvr1", &Net::nvr1, R"pbdoc(first border condition)pbdoc")
		.def_readwrite("nvr2", &Net::nvr2, R"pbdoc(second border condition)pbdoc")
		.def_readwrite("nvr3", &Net::nvr3, R"pbdoc(third border condition)pbdoc");

#ifdef VERSION_INFO
	m.attr("__version__") = VERSION_INFO;
#else
	m.attr("__version__") = "dev";
#endif
}


#endif 

/*
===============================================struct of net.txt===============================================
Nr
r1 r2 r3 ... rNr
Nz
z1 z2 z3 ...
step
hk1 kr1 hk2 kr2 ... hk(Nr-1) kr(Nr-1)
hk1 kz1 hk2 kz2 ...
Nw
W r_left z_left r_right z_right// W from 1; r_left... from 0
...
=============================================struct of border.txt==============================================
N
type NumBorder r_left z_left r_right z_right// NumBorder from 1!!!
...
==============================================struct of time.txt==============================================
Nt
t1 t2 ... tNt
stepT
hk1 kt1 hk2 kt2 ...
*/
