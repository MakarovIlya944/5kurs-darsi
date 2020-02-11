#pragma once
#include "NetGenerator.h"
#include "MatrixGenerator.h"
#include "LOS.h"

class Interface
{
	NetGenerator NetGen;
	MatrixGenerator MatGen;
	LOS Resolvent;

	char *FnetFile = "net.txt";
	char *FborderFile = "border.txt";
	char *FtimeFile = "time.txt";

	bool isMatGenReady = false;
	bool isNetGenReady = false;
	int nt, nx;
	double **Q, *ans;

	Pointd P[9]{ {4,4},{ 6,4 },{ 7,4 },{ 4,6 },{ 6,6 },{ 7,6 },{ 4,8 },{ 6,8 },{ 7,8 } };
	Pointd Pt[3]{ {4,4},{4,6},{6,6} };
	double T[3]{ 2,5,7 };
	int Nt = 3;

	double u(Pointd v, double t);

public:

	Interface()
	{
	}

	void SetNameInFiles(char *net, char *border, char *_time)
	{
		FnetFile = net;
		FborderFile = border;
		FtimeFile = _time;
	}	

	void MakeNet()
	{
		isNetGenReady = NetGen.Load(FnetFile, FborderFile, FtimeFile);
		if (isNetGenReady)
		{
			nt = NetGen.NumT;
			nx = NetGen.Num;
			Q = new double*[nt];
			ans = new double[nx];
			for (int i(0); i < nt; i++)
				Q[i] = new double[nx];
		}
	}

	void MakeMatrix()
	{
		if (isNetGenReady)
		{
			MatGen.InitMatrixGenerator(NetGen.GlobalTime, Q, NetGen.NumEl, nx);
			isMatGenReady = MatGen.GenIgJg(NetGen.nvtr);
		}
	}

	void InitFirstTime();

	void Solve()
	{
		if (isMatGenReady)
		{
			//инициализация времени
			InitFirstTime();
			for (int i(2); i < nt; i++)
			{
				MatGen.CollectGlobalMatrix(i, NetGen.nvcat, NetGen.nvtr, NetGen.GlobalNet);
				MatGen.AccountNaturalCondition(NetGen.nvr2, NetGen.nvr3, NetGen.GlobalNet);
				MatGen.AccountMainCondition(NetGen.nvr1, NetGen.GlobalNet);
				MatGen.SaveMatrix();

				Resolvent.Load();
				Resolvent.FactorLU();
				Resolvent.SolveLU();
				//Resolvent.Solve();

				ans = Resolvent.GetAnswer();
				memcpy_s(Q[i], nx * sizeof(double), ans, nx * sizeof(double));
				MatGen.ClearMatrix();
			}
		}
	}

	void Save(char *f)
	{
		FILE *F;
		fopen_s(&F, f, "w");
		for (int k(0); k < nt; k++)
		{
			fprintf_s(F, "t = %lf\n", NetGen.GlobalTime[k]);
			for (int i(0), n(NetGen.zNum); i < n; i++)
			{
				for (int j(0), nj(NetGen.rNum); j < nj; j++)
					fprintf_s(F, "%.14E;", Q[k][i*nj + j]);
				fprintf_s(F, "\n");
			}
			fprintf_s(F, "\n");
		}

		fclose(F);
	}

	void Searching()
	{
		MakeNet();
		MakeMatrix();
		Solve();
		FILE *F;
		fopen_s(&F, "ans.csv", "w");

		for (int k(0); k < nt; k++)
			fprintf_s(F, "t=%lf;8;9;10;15;16;17;22;23;24;", NetGen.GlobalTime[k]);

		fprintf_s(F, "\n");
		for (int k(0); k < nt; k++)
		{
			fprintf_s(F, "r;");
			for (int j(0); j < 9; j++)
				fprintf_s(F, "%.lf;", P[j].r);
		}

		fprintf_s(F, "\n");
		for (int k(0); k < nt; k++)
		{
			fprintf_s(F, "z;");
			for (int j(0); j < 9; j++)
				fprintf_s(F, "%.lf;", P[j].z);
		}

		fprintf_s(F, "\n");
		for (int k(0); k < nt; k++)
		{
			fprintf_s(F, "u*;");
			for (int j(0); j < 9; j++)
				fprintf_s(F, "%.lf;", u(P[j], NetGen.GlobalTime[k]));
		}

		fprintf_s(F, "\n");
		for (int k(0); k < nt; k++)
		{
			fprintf_s(F, "u(h);");
			for (int i(0), l(0), n(NetGen.zNum); i < n; i++)
				for (int j(0), nj(NetGen.rNum); j < nj; j++)
					if (abs(NetGen.GlobalNet[i*nj + j].r - P[l].r) < 1E-8 && abs(NetGen.GlobalNet[i*nj + j].z - P[l].z) < 1E-8)
					{
						fprintf_s(F, "%.14E;", Q[k][i*nj + j]);
						l++;
					}
		}
		fprintf_s(F, "\n");

		MatGen.ClearMatrix();
		SetNameInFiles("net2.txt", "border.txt","time.txt");
		MakeNet();
		MakeMatrix();
		Solve();

		for (int k(0); k < nt; k++)
		{
			fprintf_s(F, "u(h/2);");
			for (int i(0), l(0), n(NetGen.zNum); i < n; i++)
				for (int j(0), nj(NetGen.rNum); j < nj; j++)
					if (abs(NetGen.GlobalNet[i*nj + j].r - P[l].r) < 1E-8 && abs(NetGen.GlobalNet[i*nj + j].z - P[l].z) < 1E-8)
					{
						fprintf_s(F, "%.14E;", Q[k][i*nj + j]);
						l++;
					}
		}

		fprintf_s(F, "\n");

		MatGen.ClearMatrix();
		SetNameInFiles("net4.txt", "border.txt", "time.txt");
		MakeNet();
		MakeMatrix();
		Solve();

		for (int k(0); k < nt; k++)
		{
			fprintf_s(F, "u(h/4);");
			for (int i(0), l(0), n(NetGen.zNum); i < n; i++)
				for (int j(0), nj(NetGen.rNum); j < nj; j++)
					if (abs(NetGen.GlobalNet[i*nj + j].r - P[l].r) < 1E-8 && abs(NetGen.GlobalNet[i*nj + j].z - P[l].z) < 1E-8)
					{
						fprintf_s(F, "%.14E;", Q[k][i*nj + j]);
						l++;
					}
		}

		fclose(F);
	}

	void SearchingTime()
	{
		MakeNet();
		MakeMatrix();
		Solve();
		FILE *F;
		fopen_s(&F, "ansTime.csv", "w");

		fprintf_s(F, "r;4;4;4;6;6;6;7;7;7;\n");
		fprintf_s(F, "z;4;4;4;6;6;6;8;8;8;\n");

		fprintf_s(F, "t;");
		for(int i(0);i<3;i++)
			for (int j(0); j<Nt; j++)
				fprintf_s(F, "%.0lf;",T[j]);
		fprintf_s(F, "\nu*;");

		for (int j(0); j < Nt; j++)
			for (int k(0),l(0); k < nt; k++)
			if (abs(NetGen.GlobalTime[k] - T[l])<1E-8)
			{
				l++;
				fprintf_s(F, "%.lf;", u(Pt[j], NetGen.GlobalTime[k]));
			}

		fprintf_s(F, "\nu(h);");
			for (int k(0), lk(0); k < nt; k++)
			{
				if (abs(NetGen.GlobalTime[k] - T[lk]) < 1E-8)
				{
					lk++;
					fprintf_s(F, "%.14E;", Q[k][8]);
				}
			}
			for (int k(0), lk(0); k < nt; k++)
			{
				if (abs(NetGen.GlobalTime[k] - T[lk]) < 1E-8)
				{
					lk++;
					fprintf_s(F, "%.14E;", Q[k][16]);
				}
			}
			for (int k(0), lk(0); k < nt; k++)
			{
				if (abs(NetGen.GlobalTime[k] - T[lk]) < 1E-8)
				{
					lk++;
					fprintf_s(F, "%.14E;", Q[k][24]);
				}
			}

		fprintf_s(F, "\nu(h/2);");

		MatGen.ClearMatrix();
		SetNameInFiles("net.txt", "border.txt","time2.txt");
		MakeNet();
		MakeMatrix();
		Solve();

		for (int k(0), lk(0); k < nt; k++)
		{
			if (abs(NetGen.GlobalTime[k] - T[lk]) < 1E-8)
			{
				lk++;
				fprintf_s(F, "%.14E;", Q[k][8]);
			}
		}
		for (int k(0), lk(0); k < nt; k++)
		{
			if (abs(NetGen.GlobalTime[k] - T[lk]) < 1E-8)
			{
				lk++;
				fprintf_s(F, "%.14E;", Q[k][16]);
			}
		}
		for (int k(0), lk(0); k < nt; k++)
		{
			if (abs(NetGen.GlobalTime[k] - T[lk]) < 1E-8)
			{
				lk++;
				fprintf_s(F, "%.14E;", Q[k][24]);
			}
		}

		fprintf_s(F, "\nu(h/4);");

		MatGen.ClearMatrix();
		SetNameInFiles("net.txt", "border.txt", "time4.txt");
		MakeNet();
		MakeMatrix();
		Solve();

		for (int k(0), lk(0); k < nt; k++)
		{
			if (abs(NetGen.GlobalTime[k] - T[lk]) < 1E-8)
			{
				lk++;
				fprintf_s(F, "%.14E;", Q[k][8]);
			}
		}
		for (int k(0), lk(0); k < nt; k++)
		{
			if (abs(NetGen.GlobalTime[k] - T[lk]) < 1E-8)
			{
				lk++;
				fprintf_s(F, "%.14E;", Q[k][16]);
			}
		}
		for (int k(0), lk(0); k < nt; k++)
		{
			if (abs(NetGen.GlobalTime[k] - T[lk]) < 1E-8)
			{
				lk++;
				fprintf_s(F, "%.14E;", Q[k][24]);
			}
		}

		fclose(F);
	}

	double Function(double r, double z, double t)
	{
		const Pointd *p = NetGen.GlobalNet;
		const int *elems = NetGen.nvtr;
		int n = NetGen.NumEl;
		if (r < p[0].r || r > p[NetGen.Num -1].r || z < p[0].z || z > p[NetGen.Num - 1].z)
		{
			return 0;
		}

		double answer = 0;
		double R[3], Z[3], hr, hz;

		for (int i(0);i < n; i++)
		{
			if (p[elems[0]].r <= r && p[elems[0]].z <= z && p[elems[8]].r >= r && p[elems[8]].z >= z)
			{
				R[0] = (r - p[elems[1]].r)*(r - p[elems[2]].r) / (p[elems[0]].r - p[elems[1]].r) / (p[elems[0]].r - p[elems[2]].r);
				R[1] = (r - p[elems[0]].r)*(r - p[elems[2]].r) / (p[elems[1]].r - p[elems[0]].r) / (p[elems[1]].r - p[elems[2]].r);
				R[2] = (r - p[elems[0]].r)*(r - p[elems[1]].r) / (p[elems[2]].r - p[elems[0]].r) / (p[elems[2]].r - p[elems[1]].r);

				Z[0] = (z - p[elems[3]].z)*(z - p[elems[6]].z) / (p[elems[0]].z - p[elems[3]].z) / (p[elems[0]].z - p[elems[6]].z);
				Z[1] = (z - p[elems[0]].z)*(z - p[elems[6]].z) / (p[elems[3]].z - p[elems[0]].z) / (p[elems[3]].z - p[elems[6]].z);
				Z[2] = (z - p[elems[0]].z)*(z - p[elems[3]].z) / (p[elems[6]].z - p[elems[0]].z) / (p[elems[6]].z - p[elems[3]].z);
				
				for (int j(0); j < 3; j++)
					for (int k(0); k < 3; k++)
						//answer += R[k] * Z[j] * Q[t][elems[k + j * 3]];
				break;
			}
			elems += 9;
		}
		return answer;
	}

	~Interface()
	{
	}
};

