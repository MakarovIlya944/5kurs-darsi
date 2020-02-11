#pragma once
#include "NetGenerator.h"
#include <set>

using namespace std;

class MatrixGenerator
{
	bool isReady = false;
	//колличество элементов
	int nelems;
	//колличество узлов
	int npoints;
	int *ig, *jg;
	double *di, *ggl, *ggu, *f;
	const double Beta[2]{ 1,1 };
	int curType = -1;
	//Q[nt][nx]
	double **Q;
	//global time
	double *GlobalTime;
	//current time layer
	int time;

	//ѕересчитать! ст.227 соловейчик
	const double G[3][3] = {
		{ 2.33333333333333,-2.66666666666666,0.33333333333333 },
		{ -2.66666666666666,5.33333333333333,-2.66666666666666 },
		{ 0.33333333333333,-2.66666666666666,2.33333333333333 } };
	const double GH[3][3] = {
		{ 0.50000000000000,-0.66666666666666,0.16666666666666},
		{ -0.66666666666666,2.66666666666666,-2.00000000000000},
		{ 0.16666666666666,-2.00000000000000,1.83333333333333}
	};
	const double M1[3][3] = {
		{  0.1166666667,0.06666666667,-0.01666666667  },
		{ 0.06666666667,0.2666666667,0  },
		{ -0.01666666667,0,0.01666666667 } 
	};
	const double M2[3][3] = {
		{ 0.01666666667,0.,-0.01666666667 },
		{ 0.,0.2666666667,0.06666666667 },
		{ -0.01666666667,0.06666666667,0.1166666667 }
	};
	const double MH1[3][3] = {
		{ 0.01190476190,0.009523809524,-0.004761904762 },
		{ 0.009523809524,0.1142857143,0.009523809524 },
		{ -0.004761904762,0.009523809524,0.01190476190 }
	};
	const double MH2[3][3] = {
		{ 0.004761904762,-0.009523809524,-0.01190476190 },
		{ -0.009523809524,0.1523809524,0.05714285714 },
		{ -0.01190476190,0.05714285714,0.1047619048 }
	};

	set<int> *edges;

	double LocalMatrixG[9][9]{}, LocalMatrixM[9][9]{}, LocalRightPart[9];
	double LocalMatrixC[9][9];
	double LocalMatrix[9][9];

	double LocalMatrixMr[3][3];
	double LocalMatrixMz[3][3];

	int u(int i)
	{
		return i % 3 ;
	}

	int v(int i)
	{
		return i / 3 ;
	}

public:

	double u(Pointd p, double t);

	MatrixGenerator()
	{
	}

	//Constructor, that free memory for ig, init edges
	MatrixGenerator(int Nelem, int Npoints)
	{
		nelems = Nelem;
		npoints = Npoints;
		ig = new int[Npoints];
		edges = new set<int>[Npoints];
	}

	//Free memory for ig, init edges; ig - missing first 0; Q - answer
	void InitMatrixGenerator(double *&_time, double **&_Q, int Nelem, int Npoints)
	{
		nelems = Nelem;
		npoints = Npoints;
		ig = new int[Npoints];
		di = new double[Npoints]{};
		f = new double[Npoints]{};
		edges = new set<int>[Npoints];
		Q = _Q;
		GlobalTime = _time;
	}

	bool GenIgJg(int *&nvtr)
	{
		for (int i(0), N(nelems * 9); i < N; i += 9)
		{
			for (int k(0); k < 8; k++)
			{
				int tmp = nvtr[i + k];
				for (int j(k + 1); j < 9; j++)
					edges[nvtr[i + j]].insert(tmp);
			}
		}
		int nJg = 0;
		ig[0] = 0;
		for (int i(1); i < npoints; i++)
		{
			ig[i] = ig[i - 1] + edges[i].size();
			nJg += ig[i];
		}
		jg = new int[nJg];
		for (int i(0), k(0); i < npoints; i++)
			for each (int var in edges[i])
				jg[k++] = var;
		ggl = new double[nJg] {};
		ggu = new double[nJg] {};
		isReady = true;
		return true;
	}

	double CoefficientDiffusion(Pointd v);

	double CoefficientSigma(Pointd v);

	void MakeLocalMatrixMassR(double r, double g1, double g2, double hr)
	{
		for (int i(0); i < 3; i++)
			for (int j(0); j < 3; j++)
				LocalMatrixMr[i][j] = hr*(r * (M1[i][j] * g1 + M2[i][j] * g2) + hr * (MH1[i][j] * g1 + MH2[i][j] * g2));
	}

	void MakeLocalMatrixMassZ(double g1, double g2, double hz)
	{
		for (int i(0); i < 3; i++)
			for (int j(0); j < 3; j++)
				LocalMatrixMz[i][j] = hz * (M1[i][j] * g1 + M2[i][j] * g2);
	}

	void MakeLocalMatrix(int *&elements, Pointd *&globnet)
	{
		double hr = globnet[elements[2]].r - globnet[elements[0]].r;
		double hz = globnet[elements[6]].z - globnet[elements[0]].z;
		double lambda = CoefficientDiffusion(globnet[elements[4]]), r = globnet[elements[0]].r;
		double dt0 = GlobalTime[time] - GlobalTime[time - 1], dt1 = GlobalTime[time - 1] - GlobalTime[time - 2], dt2 = GlobalTime[time] - GlobalTime[time - 2];
		double k1 = (dt2 + dt0) / (dt2 * dt0), k2 = dt0 / (dt2 * dt1), k3 = dt2 / (dt1 * dt0);

		MakeLocalMatrixMassR(r, CoefficientSigma(globnet[elements[3]]), CoefficientSigma(globnet[elements[5]]), hr);
		MakeLocalMatrixMassZ(CoefficientSigma(globnet[elements[1]]), CoefficientSigma(globnet[elements[7]]), hz);
		int vi, vj, ui, uj;
		for (int i(0); i < 9; i++)
		{
			vi = v(i); ui = u(i);
			for (int j(0); j < 9; j++)
			{
				vj = v(j); uj = u(j);
				LocalMatrixG[i][j] = (r * G[ui][uj] / hr + GH[ui][uj]) * (hz * (M1[vi][vj] + M2[vi][vj]));
				LocalMatrixG[i][j] += G[vi][vj] / hz * hr*(r * (M1[ui][uj] + M2[ui][uj]) + hr * (MH1[ui][uj] + MH2[ui][uj]));
				LocalMatrixM[i][j] = LocalMatrixMr[ui][uj] * LocalMatrixMz[vi][vj];
				LocalMatrix[i][j] = lambda * LocalMatrixG[i][j] + LocalMatrixM[i][j] * k1;

				f[elements[i]] += RightPart(globnet[elements[j]], GlobalTime[time]) * (((M1[ui][uj] + M2[ui][uj]) * r + (MH1[ui][uj] + MH2[ui][uj])*hr) *(M1[vi][vj] + M2[vi][vj])* hr*hz);
				f[elements[i]] += LocalMatrixM[i][j]*(k3*Q[time - 1][elements[j]] - k2*Q[time - 2][elements[j]]);
			}
		}
	}

	//elements - nvtr, GlobalPoints
	bool CollectGlobalMatrix(int _t, int *&materials, int *&elements, Pointd *&globnet)
	{
		time = _t;
		if (!isReady)
			return false;
		int *pointelem = elements;
		for (int elem(0); elem < nelems; elem++)
		{
			curType = materials[elem];
			MakeLocalMatrix(pointelem, globnet);
			//ст.500 Cоловейчик 
			//add diag elements
			for (int i(0); i < 9; i++)
				di[pointelem[i]] += LocalMatrix[i][i];

			for (int i(0), ibeg; i < 9; i++)
			{
				ibeg = ig[pointelem[i] < 1 ? 0 : pointelem[i] - 1];
				for (int j(0); j < i; j++)
				{
					int iend = ig[pointelem[i]] - 1;
					while (jg[ibeg] != pointelem[j])//binary find
					{
						int ind = (ibeg + iend) % 2 ? (ibeg + iend) / 2 + 1 : (ibeg + iend) / 2;
						(jg[ind] <= pointelem[j] ? ibeg : iend) = ind; //if (jg[ind] <= pointelem[j])   ibeg = ind;   else   iend = ind;
					}
					ggl[ibeg] += LocalMatrix[i][j];
					ggu[ibeg] += LocalMatrix[j][i];
					ibeg++;
				}
			}
			pointelem += 9;
		}
		return true;
	}

	//Дописать
	double RightPart(Pointd v,double _t);

	double FirstCond(Pointd v);

	double SecondCond(Pointd v);

	double ThirdCond(Pointd v);

	void AccountNaturalCondition(list<int> &nvr2, list<int> &nvr3, Pointd *&globnet)
	{
		double h,SecCond_a, SecCond_b, SecCond_c;
		curType = -2;
		int var;
		if (!nvr2.empty())
			for (auto c = nvr2.begin(), a = c++, b = c++;;)
			{
				var = *a;
				if (var == -1)
				{
					curType = -2;
					a++;
					b++;
					c++;
				}
				else if (curType == -2)
				{
					curType = var;
					a++;
					b++;
					c++;
					if ((*c) - (*a) == 2)//r
						h = globnet[(*c)].r - globnet[(*a)].r;
					else//z
						h = globnet[(*c)].z - globnet[(*a)].z;
				}
				else
				{
					SecCond_a = SecondCond(globnet[*a]);
					SecCond_b = SecondCond(globnet[*b]);
					SecCond_c = SecondCond(globnet[*c]);
					if ((*c) - (*a) == 2)//r
					{
						MakeLocalMatrixMassR(globnet[*a].r, CoefficientSigma(globnet[*a]), CoefficientSigma(globnet[*c]), h);
						f[*a] += SecCond_a * LocalMatrixMr[0][0] + SecCond_b * LocalMatrixMr[0][1] + SecCond_c * LocalMatrixMr[0][2];
						f[*b] += SecCond_a * LocalMatrixMr[1][0] + SecCond_b * LocalMatrixMr[1][1] + SecCond_c * LocalMatrixMr[1][2];
						f[*c] += SecCond_a * LocalMatrixMr[2][0] + SecCond_b * LocalMatrixMr[2][1] + SecCond_c * LocalMatrixMr[2][2];
					}
					else//z
					{
						MakeLocalMatrixMassZ(CoefficientSigma(globnet[*a]), CoefficientSigma(globnet[*c]), h);
						f[*a] += SecCond_a * LocalMatrixMz[0][0] + SecCond_b * LocalMatrixMz[0][1] + SecCond_c * LocalMatrixMz[0][2];
						f[*b] += SecCond_a * LocalMatrixMz[1][0] + SecCond_b * LocalMatrixMz[1][1] + SecCond_c * LocalMatrixMz[1][2];
						f[*c] += SecCond_a * LocalMatrixMz[2][0] + SecCond_b * LocalMatrixMz[2][1] + SecCond_c * LocalMatrixMz[2][2];
					}
					a++;
					a++;
					b++;
					b++;
					if (*b == -1)
					{
						a++;
						b++;
						c++;
					}
					c++;
					if (c == nvr2.end())
						break;
					c++;
				}
			}

		curType = -2;

		double koef = 0;
		int indexis[3];
		if (!nvr3.empty())
			for (auto c = nvr3.begin(), a = c++, b = c++;;)
			{
				indexis[0] = *a;
				indexis[1] = *b;
				indexis[2] = *c;
				if (indexis[0] == -1)
				{
					curType = -2;
					a++;
					b++;
					c++;
				}
				else if (curType == -2)
				{
					curType = indexis[0];
					a++;
					b++;
					c++;
					if ((*c) - indexis[0] == 2)//r
						h = globnet[indexis[2]].r - globnet[indexis[0]].r;
					else//z
						h = globnet[indexis[2]].z - globnet[indexis[0]].z;
				}
				else
				{
					SecCond_a = SecondCond(globnet[*a]);
					SecCond_b = SecondCond(globnet[*b]);
					SecCond_c = SecondCond(globnet[*c]);
					if ((*c) - (*a) == 2)//r
					{
						MakeLocalMatrixMassR(globnet[*a].r, CoefficientSigma(globnet[*a]), CoefficientSigma(globnet[*c]), h);
						f[*a] += Beta[curType] * SecCond_a * LocalMatrixMr[0][0] + SecCond_b * LocalMatrixMr[0][1] + SecCond_c * LocalMatrixMr[0][2];
						f[*b] += Beta[curType] * SecCond_a * LocalMatrixMr[1][0] + SecCond_b * LocalMatrixMr[1][1] + SecCond_c * LocalMatrixMr[1][2];
						f[*c] += Beta[curType] * SecCond_a * LocalMatrixMr[2][0] + SecCond_b * LocalMatrixMr[2][1] + SecCond_c * LocalMatrixMr[2][2];

						for (int i(0); i < 3; i++)
							di[indexis[i]] += Beta[curType] * LocalMatrixMr[i][i];

						for (int i(0); i < 3; i++)//local matrix
						{
							int ibeg = ig[indexis[i] < 1 ? 0 : indexis[i] - 1];
							for (int j(0); j < i; j++)
							{
								int iend = ig[indexis[i]] + 1;
								while (jg[ibeg] != indexis[j])//binary find
								{
									int ind = (ibeg + iend) % 2 ? (ibeg + iend) / 2 + 1 : (ibeg + iend) / 2;
									if (jg[ind] <= indexis[j])
										ibeg = ind;
									else
										iend = ind;
								}
								ggl[ibeg] += koef * LocalMatrixMr[i][j];
								ggu[ibeg] += koef * LocalMatrixMr[j][i];
								ibeg++;
							}
						}
					}
					else//z
					{
						MakeLocalMatrixMassZ(CoefficientSigma(globnet[*a]), CoefficientSigma(globnet[*c]), h);
						f[*a] += Beta[curType] * SecCond_a * LocalMatrixMz[0][0] + SecCond_b * LocalMatrixMz[0][1] + SecCond_c * LocalMatrixMz[0][2];
						f[*b] += Beta[curType] * SecCond_a * LocalMatrixMz[1][0] + SecCond_b * LocalMatrixMz[1][1] + SecCond_c * LocalMatrixMz[1][2];
						f[*c] += Beta[curType] * SecCond_a * LocalMatrixMz[2][0] + SecCond_b * LocalMatrixMz[2][1] + SecCond_c * LocalMatrixMz[2][2];

						for (int i(0); i < 3; i++)
							di[indexis[i]] += Beta[curType] * LocalMatrixMz[i][i];

						for (int i(0); i < 3; i++)//local matrix
						{
							int ibeg = ig[indexis[i] < 1 ? 0 : indexis[i] - 1];
							for (int j(0); j < i; j++)
							{
								int iend = ig[indexis[i]] + 1;
								while (jg[ibeg] != indexis[j])//binary find
								{
									int ind = (ibeg + iend) % 2 ? (ibeg + iend) / 2 + 1 : (ibeg + iend) / 2;
									if (jg[ind] <= indexis[j])
										ibeg = ind;
									else
										iend = ind;
								}
								ggl[ibeg] += koef * LocalMatrixMz[i][j];
								ggu[ibeg] += koef * LocalMatrixMz[j][i];
								ibeg++;
							}
						}
					}

					a++; a++;
					b++; b++;
					if (*b == -1)
					{
						a++;
						b++;
						c++;
					}
					c++;
					if (c == nvr3.end())
						break;
					c++;
				}
			}
	}

	void AccountMainCondition(list<int> &nvr1, Pointd *&globnet)
	{
		curType = 0;
		int  ibeg, iend;
		bool flag = true;
		double cond;
		for each (int var in nvr1)
		{
			if (var == -1)
				curType = 0;
			else if (curType == 0)
				curType = var;
			else
			{
				cond = FirstCond(globnet[var]);
				di[var] = 1;
				f[var] = cond;
				ibeg = var ? ig[var - 1] : 0;
				iend = ig[var];
				for (int i(ibeg), n(iend); i < n; i++)
				{
					ggl[i] = 0;
					//f[jg[i]] -= cond * ggu[i];
					//ggu[i] = 0;
				}

				for (int i(var + 1), n(npoints), tmp=-1, ind; i < n; i++)
				{
					iend = ig[i];
					ibeg = ig[i - 1];
					if (jg[iend-1] < var || jg[ibeg] > var)
						continue;
					while (jg[ibeg] != var && flag)//binary find
					{
						ind = (ibeg + iend) % 2 ? (ibeg + iend) / 2 + 1 : (ibeg + iend) / 2;
						flag = ind != tmp;
						if (jg[ind] <= var)
							ibeg = ind;
						else
							iend = ind;
						tmp = ind;
					}
					if (flag)//если нашли индекс в ibeg
					{
						ggu[ibeg] = 0;
						//f[GetI(ibeg)+1] -= ggl[ibeg] * cond;
						//ggl[ibeg] = 0;
					}
					else
						flag = true;
				}
			}
		}
	}

	void SaveMatrix()
	{
		ofstream out1, out2, out3, out4;

		out1.precision(15);
		out2.precision(15);
		out3.precision(15);
		out4.precision(15);

		out1.open("kuslau.txt");
		out1 << npoints << endl << 10000 << endl << 1E-14;
		out1.close();

		out1.open("x.txt");//начальное
		out2.open("pr.txt");
		out3.open("di.txt");
		out4.open("ig.txt");
		out4 << 0 << ' ';
		for (int i(0); i < npoints; i++)
		{
			out1 << 0 << ' ';
			out2 << f[i] << ' ';
			out3 << di[i] << ' ';
			out4 << ig[i] << ' ';
		}
		out1.close();
		out2.close();
		out3.close();
		out4.close();

		int n = ig[npoints - 1];

		out1.open("jg.txt");//начальное
		out2.open("ggl.txt");
		out3.open("ggu.txt");
		for (int i(0); i < n; i++)
		{
			out1 << jg[i] << ' ';
			out2 << ggl[i] << ' ';
			out3 << ggu[i] << ' ';
		}
		out1.close();
		out2.close();
		out3.close();
	}

	void ClearMatrix()
	{
		std::fill(ggl, ggl + ig[npoints-1], 0);
		std::fill(ggu, ggu + ig[npoints-1], 0);
		std::fill(di, di + npoints, 0);
		std::fill(f, f + npoints, 0);
	}

	~MatrixGenerator()
	{
		delete ig;
		delete jg;
		delete ggl;
		delete ggu;
		delete di;
		delete f;
	}
};