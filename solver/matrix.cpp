#include "matrix.hpp"
int curType = 0;
double CoefficientDiffusion(Pointd v)
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

double CoefficientSigma(Pointd v)
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

double RightPart(Pointd v, double t)
{
	return 6 * t * t;
}

double u(Pointd v, double t)
{
	return t * t * t;
}

double FirstCond(Pointd v, double t)
{
	return u(v, t);
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
		return v.z * v.z;
		break;
	case 3:
		//return v.r*v.r + 25;
		//return v.r*v.r*v.r+6*v.r+632;
		return 25;
		break;
	case 4:
		//return v.z*v.z + 1;
		//return v.z*v.z*v.z*v.z+2*v.z+4;
		return v.z * v.z;
		break;
	default:
		return 0;
		break;
	}
}

double SecondCond(Pointd v, double t)
{
	return -2 * t * t;
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
		return 2 * (4 * v.z * v.z * v.z + v.r + 1);
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

double ThirdCond(Pointd v)
{
	switch (curType)
	{
	case 1:
		//return -2*v.z;
		return v.r * v.r + 1 - v.z;
		break;
	case 2:
		//return 2*v.r;
		return 2 * (3 * v.r * v.r + v.z + 1);
		break;
	case 3:
		//return 2*v.z;
		return 2 * (4 * v.z * v.z * v.z + v.r + 1);
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


MatrixGenerator::MatrixGenerator(Net *_net, double **&_Q)
{
	net = _net;
	edges = new set<int>[net->Num];
	ig = new int[net->Num];
	di = new double[net->Num]{};
	f = new double[net->Num]{};
	matricies = map<double, Matrix>();
	Q = _Q;

	for (int i(0), N(net->NumEl * 9); i < N; i += 9)
	{
		for (int k(0); k < 8; k++)
		{
			int tmp = net->nvtr[i + k];
			for (int j(k + 1); j < 9; j++)
				edges[net->nvtr[i + j]].insert(tmp);
		}
	}
	int nJg = 0;
	ig[0] = 0;
	for (int i(1); i < net->Num; i++)
	{
		ig[i] = ig[i - 1] + edges[i].size();
		nJg += ig[i];
	}
	jg = new int[nJg];
	for (int i(0), k(0); i < net->Num; i++)
		for (int var : edges[i])
			jg[k++] = var;
	ggl = new double[nJg]{};
	ggu = new double[nJg]{};
}

void MatrixGenerator::MakeLocalMatrixMassR(double r, double g1, double g2, double hr)
{
	for (int i(0); i < 3; i++)
		for (int j(0); j < 3; j++)
			LocalMatrixMr[i][j] = hr * (r * (M1[i][j] * g1 + M2[i][j] * g2) + hr * (MH1[i][j] * g1 + MH2[i][j] * g2));
}

void MatrixGenerator::MakeLocalMatrixMassZ(double g1, double g2, double hz)
{
	for (int i(0); i < 3; i++)
		for (int j(0); j < 3; j++)
			LocalMatrixMz[i][j] = hz * (M1[i][j] * g1 + M2[i][j] * g2);
}

void MatrixGenerator::MakeLocalMatrix(int time, vector<int>::iterator elements)
{
	vector<Pointd> &globnet = net->GlobalNet;
	vector<double> &globtime = net->GlobalTime;
	double hr = globnet[*(elements + 2)].r - globnet[*elements].r;
	double hz = globnet[*(elements + 6)].z - globnet[*elements].z;
	double lambda = CoefficientDiffusion(globnet[elements[4]]), r = globnet[elements[0]].r;
	double dt0 = globtime[time] - globtime[time - 1], dt1 = globtime[time - 1] - globtime[time - 2], dt2 = globtime[time] - globtime[time - 2];
	double k1 = (dt2 + dt0) / (dt2 * dt0), k2 = dt0 / (dt2 * dt1), k3 = dt2 / (dt1 * dt0);

	MakeLocalMatrixMassR(r, CoefficientSigma(globnet[elements[3]]), CoefficientSigma(globnet[elements[5]]), hr);
	MakeLocalMatrixMassZ(CoefficientSigma(globnet[elements[1]]), CoefficientSigma(globnet[elements[7]]), hz);
	int vi, vj, ui, uj;
	for (int i(0); i < 9; i++)
	{
		vi = v(i);
		ui = u(i);
		for (int j(0); j < 9; j++)
		{
			vj = v(j);
			uj = u(j);
			LocalMatrixG[i][j] = (r * G[ui][uj] / hr + GH[ui][uj]) * (hz * (M1[vi][vj] + M2[vi][vj]));
			LocalMatrixG[i][j] += G[vi][vj] / hz * hr * (r * (M1[ui][uj] + M2[ui][uj]) + hr * (MH1[ui][uj] + MH2[ui][uj]));
			LocalMatrixM[i][j] = LocalMatrixMr[ui][uj] * LocalMatrixMz[vi][vj];
			LocalMatrix[i][j] = lambda * LocalMatrixG[i][j] + LocalMatrixM[i][j] * k1;

			f[elements[i]] += RightPart(globnet[elements[j]], globtime[time]) * (((M1[ui][uj] + M2[ui][uj]) * r + (MH1[ui][uj] + MH2[ui][uj]) * hr) * (M1[vi][vj] + M2[vi][vj]) * hr * hz);
			f[elements[i]] += LocalMatrixM[i][j] * (k3 * Q[time - 1][elements[j]] - k2 * Q[time - 2][elements[j]]);
		}
	}
}

void MatrixGenerator::CollectGlobalMatrix()
{
	for (int time(0), ntime(net->GlobalTime.size()); time < ntime; time++)
	{
		auto pointelem = net->nvtr.begin();
		for (int elem(0); elem < net->NumEl; elem++)
		{
			curType = net->nvcat[elem];
			MakeLocalMatrix(time, pointelem);
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
					while (jg[ibeg] != pointelem[j]) //binary find
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
		AccounConditions(time);
	}
}

void MatrixGenerator::AccounConditions(int time)
{
	AccountNaturalCondition(time);
	AccountMainCondition(time);
}

void MatrixGenerator::AccountNaturalCondition(int time)
{
	double h, SecCond_a, SecCond_b, SecCond_c;
	curType = -2;
	vector<Pointd> &globnet = net->GlobalNet;
	int var;
	if (net->nvr2.empty())
		for (auto c = net->nvr2.begin(), a = c++, b = c++;;)
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
				if ((*c) - (*a) == 2) //r
					h = globnet[(*c)].r - globnet[(*a)].r;
				else //z
					h = globnet[(*c)].z - globnet[(*a)].z;
			}
			else
			{
				SecCond_a = SecondCond(globnet[*a], time);
				SecCond_b = SecondCond(globnet[*b], time);
				SecCond_c = SecondCond(globnet[*c], time);
				if ((*c) - (*a) == 2) //r
				{
					MakeLocalMatrixMassR(globnet[*a].r, CoefficientSigma(globnet[*a]), CoefficientSigma(globnet[*c]), h);
					f[*a] += SecCond_a * LocalMatrixMr[0][0] + SecCond_b * LocalMatrixMr[0][1] + SecCond_c * LocalMatrixMr[0][2];
					f[*b] += SecCond_a * LocalMatrixMr[1][0] + SecCond_b * LocalMatrixMr[1][1] + SecCond_c * LocalMatrixMr[1][2];
					f[*c] += SecCond_a * LocalMatrixMr[2][0] + SecCond_b * LocalMatrixMr[2][1] + SecCond_c * LocalMatrixMr[2][2];
				}
				else //z
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
				if (c == net->nvr2.end())
					break;
				c++;
			}
		}

	curType = -2;

	double koef = 0;
	int indexis[3];
	if (!(net->nvr3.empty()))
		for (auto c = net->nvr3.begin(), a = c++, b = c++;;)
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
				if ((*c) - indexis[0] == 2) //r
					h = globnet[indexis[2]].r - globnet[indexis[0]].r;
				else //z
					h = globnet[indexis[2]].z - globnet[indexis[0]].z;
			}
			else
			{
				SecCond_a = SecondCond(globnet[*a], time);
				SecCond_b = SecondCond(globnet[*b], time);
				SecCond_c = SecondCond(globnet[*c], time);
				if ((*c) - (*a) == 2) //r
				{
					MakeLocalMatrixMassR(globnet[*a].r, CoefficientSigma(globnet[*a]), CoefficientSigma(globnet[*c]), h);
					f[*a] += Beta[curType] * SecCond_a * LocalMatrixMr[0][0] + SecCond_b * LocalMatrixMr[0][1] + SecCond_c * LocalMatrixMr[0][2];
					f[*b] += Beta[curType] * SecCond_a * LocalMatrixMr[1][0] + SecCond_b * LocalMatrixMr[1][1] + SecCond_c * LocalMatrixMr[1][2];
					f[*c] += Beta[curType] * SecCond_a * LocalMatrixMr[2][0] + SecCond_b * LocalMatrixMr[2][1] + SecCond_c * LocalMatrixMr[2][2];

					for (int i(0); i < 3; i++)
						di[indexis[i]] += Beta[curType] * LocalMatrixMr[i][i];

					for (int i(0); i < 3; i++) //local matrix
					{
						int ibeg = ig[indexis[i] < 1 ? 0 : indexis[i] - 1];
						for (int j(0); j < i; j++)
						{
							int iend = ig[indexis[i]] + 1;
							while (jg[ibeg] != indexis[j]) //binary find
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
				else //z
				{
					MakeLocalMatrixMassZ(CoefficientSigma(globnet[*a]), CoefficientSigma(globnet[*c]), h);
					f[*a] += Beta[curType] * SecCond_a * LocalMatrixMz[0][0] + SecCond_b * LocalMatrixMz[0][1] + SecCond_c * LocalMatrixMz[0][2];
					f[*b] += Beta[curType] * SecCond_a * LocalMatrixMz[1][0] + SecCond_b * LocalMatrixMz[1][1] + SecCond_c * LocalMatrixMz[1][2];
					f[*c] += Beta[curType] * SecCond_a * LocalMatrixMz[2][0] + SecCond_b * LocalMatrixMz[2][1] + SecCond_c * LocalMatrixMz[2][2];

					for (int i(0); i < 3; i++)
						di[indexis[i]] += Beta[curType] * LocalMatrixMz[i][i];

					for (int i(0); i < 3; i++) //local matrix
					{
						int ibeg = ig[indexis[i] < 1 ? 0 : indexis[i] - 1];
						for (int j(0); j < i; j++)
						{
							int iend = ig[indexis[i]] + 1;
							while (jg[ibeg] != indexis[j]) //binary find
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
				if (c == net->nvr3.end())
					break;
				c++;
			}
		}
}

void MatrixGenerator::AccountMainCondition(int time)
{
	curType = 0;
	int ibeg, iend;
	bool flag = true;
	double cond;
	for (int var : net->nvr1)
	{
		if (var == -1)
			curType = 0;
		else if (curType == 0)
			curType = var;
		else
		{
			cond = FirstCond(net->GlobalNet[var], time);
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

			for (int i(var + 1), n(net->Num), tmp = -1, ind; i < n; i++)
			{
				iend = ig[i];
				ibeg = ig[i - 1];
				if (jg[iend - 1] < var || jg[ibeg] > var)
					continue;
				while (jg[ibeg] != var && flag) //binary find
				{
					ind = (ibeg + iend) % 2 ? (ibeg + iend) / 2 + 1 : (ibeg + iend) / 2;
					flag = ind != tmp;
					if (jg[ind] <= var)
						ibeg = ind;
					else
						iend = ind;
					tmp = ind;
				}
				if (flag) //если нашли индекс в ibeg
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

void MatrixGenerator::SaveMatrix()
{
	ofstream out1, out2, out3, out4;

	out1.precision(15);
	out2.precision(15);
	out3.precision(15);
	out4.precision(15);

	out1.open("kuslau.txt");
	out1 << net->Num << endl
		 << 10000 << endl
		 << 1E-14;
	out1.close();

	out1.open("x.txt"); //начальное
	out2.open("pr.txt");
	out3.open("di.txt");
	out4.open("ig.txt");
	out4 << 0 << ' ';
	for (int i(0); i < net->Num; i++)
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

	int n = ig[net->Num - 1];

	out1.open("jg.txt"); //начальное
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

void MatrixGenerator::ClearMatrix()
{
	std::fill(ggl, ggl + ig[net->Num - 1], 0);
	std::fill(ggu, ggu + ig[net->Num - 1], 0);
	std::fill(di, di + net->Num, 0);
	std::fill(f, f + net->Num, 0);
}

MatrixGenerator::~MatrixGenerator()
{
	delete ig;
	delete jg;
	delete ggl;
	delete ggu;
	delete di;
	delete f;
}