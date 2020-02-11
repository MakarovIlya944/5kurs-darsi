#pragma once
#include <fstream>
#include <iostream>
#include <list>
#include <iterator>

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

class NetGenerator
{
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
	Pointd *GlobalNet;
	//time net
	double *GlobalTime;
	//nvtr - elements
	int *nvtr;
	//nvcat - materials
	int *nvcat;
	//nvr1 - first border condition
	list<int> nvr1;
	//nvr2 - second border condition
	list<int> nvr2;
	//nvr3 - third border condition
	list<int> nvr3;
	//isReady - true if data ready
	bool isReady = false;

	NetGenerator()
	{

	}

	//FNameN - net.txt , FNameB - border.txt , FNameT - time.txt
	bool Load(const char*FNameN, const char*FNameB, const char*FNameT)
	{
		ifstream in(FNameN);

		int nr, nz;
		in >> nr;
		double *R = new double[nr];
		for (int i(0); i < nr; i++)
			in >> R[i];

		in >> nz;
		double *Z = new double[nz];
		for (int i(0); i < nz; i++)
			in >> Z[i];

		double step;
		in >> step;
		step *= 0.5;
		double *kr = new double[nr - 1], *kz = new double[nz - 1];
		int *hr = new int[nr - 1], nhr(0), *hz = new int[nz - 1], nhz(0);
		for (int i(0); i < nr - 1; i++)
		{
			in >> hr[i];
			in >> kr[i];
			hr[i] *= 2;
			nhr += hr[i];
		}
		for (int i(0); i < nz - 1; i++)
		{
			in >> hz[i];
			in >> kz[i];
			hz[i] *= 2;
			nhz += hz[i];
		}

		int nw;
		in >> nw;
		int *W = new int[nw];
		Pointi *area = new Pointi[nw * 2]{};
		for (int i(0), j(0); i < nw; i++)
		{
			in >> W[i];
			in >> area[j].r;
			in >> area[j++].z;
			in >> area[j].r;
			in >> area[j++].z;
		}

		double *Rh = new double[nhr + 1], *Zh = new double[nhz + 1];
		for (int i(0), k(0); i < nr - 1; i++)
		{
			for (int j(0); j < hr[i]; j++, k++)
				Rh[k] = R[i] + j*step*kr[i];
			hr[i] += i != 0 ? hr[i - 1] : 0;
		}
		Rh[nhr] = R[nr - 1];

		for (int i(0), k(0); i < nz - 1; i++)
		{
			for (int j(0); j < hz[i]; j++, k++)
				Zh[k] = Z[i] + j*step*kz[i];
			hz[i] += i != 0 ? hz[i - 1] : 0;
		}
		Zh[nhz] = Z[nz - 1];

		delete R;
		delete Z;

		rNum = nhr + 1;
		zNum = nhz + 1;
		Num = rNum*zNum;
		GlobalNet = new Pointd[Num];
		int n = nhr*nhz / 4;
		nvtr = new int[n * 9];
		nvcat = new int[n];
		NumEl = n;
		double tmp;
		for (int i(0), w(0), currentelem(0); i < zNum; i++)
		{
			tmp = Zh[i];
			for (int j(0); j < rNum; j++)
			{
				if ((i % 2) && (j % 2))
				{
					for (int k(0), n(nw * 2); k < n && w == 0; k += 2)
					{
						if ((area[k].r == 0 || hr[area[k].r - 1] < j) && j < hr[area[k + 1].r - 1])
							if (i < hz[area[k + 1].z - 1] && (area[k].z == 0 || hz[area[k].z - 1] < i))
								w = W[k / 2];
					}
					nvcat[currentelem++] = w;
					w = 0;
				}
				GlobalNet[i*rNum + j] = Pointd(Rh[j], tmp);
			}
		}

		GenerationNVTR();
		in.close();

		in = ifstream(FNameB);

		nvr1.clear();
		nvr2.clear();
		nvr3.clear();

		int nb;
		in >> nb;
		int *B = new int[nb], *V = new int[nb];
		Pointi *border = new Pointi[nb * 2]{};
		for (int i(0), j(0); i < nb; i++)
		{
			in >> B[i];
			in >> V[i];
			in >> border[j].r;
			in >> border[j++].z;
			in >> border[j].r;
			in >> border[j++].z;
		}

		for (int i(0), j(0); i < nb; i++, j++)
		{
			list<int> *iterator = &nvr1;

			switch (B[i])
			{
			case 2:
				iterator = &nvr2;
				break;
			case 3:
				iterator = &nvr3;
			}
			iterator->push_back(V[i]);
			if (border[j].r == border[j + 1].r)
			{
				int l = border[j].z != 0 ? hz[border[j].z - 1] : 0;
				int k = border[j].r != 0 ? hr[border[j].r - 1] : 0;
				if(border[j].z == border[j + 1].z)
					iterator->push_back(l*rNum + k);
				int nl(hz[border[++j].z - 1]);
				for (; l <= nl; l++)
					iterator->push_back(l*rNum + k);
			}
			else
			{
				int l = border[j].r != 0 ? hr[border[j].r - 1] : 0;
				int k = border[j].z != 0 ? hz[border[j].z - 1] : 0;
				int nl(hr[border[++j].r - 1]);
				for (; l <= nl; l++)
					iterator->push_back(k*rNum + l);
			}
			iterator->push_back(-1);
		}

		in.close();
		in.open(FNameT);

		int nt;
		in >> nt;
		double *T = new double[nt];
		for (int i(0); i < nt; i++)
			in >> T[i];

		in >> step;

		double *kt = new double[nt - 1];
		int *ht = new int[nt - 1], nht(0);
		for (int i(0); i < nt - 1; i++)
		{
			in >> ht[i];
			in >> kt[i];
			nht += ht[i];
		}

		GlobalTime = new double[nht + 1];
		GlobalTime[0] = T[0];
		for (int i(0),k(0); k < nht; i++)
		{
			tmp = step * kt[i];
			for (int j(0); j < ht[i]; j++,k++)
				GlobalTime[k + 1] = GlobalTime[k] + tmp;
		}
		NumT = nht + 1;
		in.close();


		delete kt;
		delete ht;
		delete T;
		delete Rh;
		delete Zh;
		delete hz;
		delete hr;
		delete kr;
		delete kz;
		delete W;
		delete border;
		
		isReady = true;
		return isReady;
	}

	void GenerationNVTR()
	{
		int zn = zNum - 1, rn = rNum - 1;
		for (int z(0), current(0); z < zn; z += 2)
			for (int r(0); r < rn; r += 2)
			{
				int index = r + z*rNum;
				nvtr[current++] = index;
				nvtr[current++] = index + 1;
				nvtr[current++] = index + 2;

				nvtr[current++] = index + rNum;
				nvtr[current++] = index + rNum + 1;
				nvtr[current++] = index + rNum + 2;

				nvtr[current++] = index + 2 * rNum;
				nvtr[current++] = index + 2 * rNum + 1;
				nvtr[current++] = index + 2 * rNum + 2;
			}
	}

	void SaveNet(char *file)
	{
		ofstream out;
		out.open(file);

		for (int i(0); i < NumT; i++)
			out << GlobalTime[i];
		for (int i(0); i < rNum; i++)
			out << GlobalNet[i].r;
		for (int i(0), j(0); i < zNum; i++, j += rNum)
			out << GlobalNet[i].z;


		out.close();
	}

	~NetGenerator()
	{
		delete GlobalNet;
		delete nvcat;
		delete nvtr;
		nvr1.clear();
		nvr2.clear();
		nvr3.clear();
	}
};

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
