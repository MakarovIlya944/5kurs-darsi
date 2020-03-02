//#pragma once 
#include <iostream>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include<conio.h>
#include<set>
#include<list>

using namespace std;

//ThreeStepMethod
class LOS
{
	double normar;

	int n, maxiter, iter;
	double eps;

	double *gglLU, *gguLU, *diLU;
	double *gglA, *gguA, *diA;
	double *diD;
	int *ig, *jg;

	double *r, *z, *x, *p, *pr, *Ax, *Ay, *x0;

public:
	//c=b-a
	void sub(double *a, double *b, double *c)
	{
		for (int i(0); i < n; i++)
			c[i] = b[i] - a[i];
	}

	//return (a,b)
	double scalar(double *a, double *b)
	{
		double ans = 0;
		for (int i(0); i < n; i++)
			ans += a[i] * b[i];
		return ans;
	}

	//d=a+b*c
	void addmult(double *a, double b, double *c, double *d)
	{
		for (int i(0); i < n; i++)
			d[i] = a[i] + b*c[i];
	}

	double* GetAnswer()
	{
		return x;
	}

	void Load()
	{
		ifstream in;
		in.open("kuslau.txt");
		in >> n;
		in >> maxiter;
		in >> eps;
		in.close();

		diLU = new double[n] {};
		diA = new double[n] {};
		diD = new double[n] {};
		ig = new int[n + 1]{};
		pr = new double[n] {};
		r = new double[n] {};
		z = new double[n] {};
		x = new double[n] {};
		p = new double[n] {};
		Ax = new double[n] {};
		Ay = new double[n] {};
		x0 = new double[n] {};

		in.open("x.txt");
		for (int i(0); i < n; i++)
		{
			in >> x[i];
			x0[i] = x[i];
		}
		in.close();

		in.open("pr.txt");
		for (int i(0); i < n; i++)
			in >> pr[i];
		in.close();

		in.open("di.txt");
		for (int i(0); i < n; i++)
			in >> diA[i];
		in.close();

		in.open("ig.txt");
		for (int i(0); i <= n; i++)
		{
			in >> ig[i];
			//ig[i]--;
		}
		in.close();

		int tmp = ig[n];
		jg = new int[tmp];
		gglA = new double[tmp] {};
		gguA = new double[tmp] {};
		gglLU = new double[tmp] {};
		gguLU = new double[tmp] {};

		in.open("jg.txt");
		for (int i(0); i < tmp; i++)
		{
			in >> jg[i];
			//jg[i]--;
		}
		in.close();

		in.open("ggl.txt");
		for (int i(0); i < tmp; i++)
			in >> gglA[i];
		in.close();

		in.open("ggu.txt");
		for (int i(0); i < tmp; i++)
			in >> gguA[i];
		in.close();
	}

	void Save(char *file)
	{
		FILE *f;
		fopen_s(&f, file, "w");
		for (int i(0); i < n; i++)
			fprintf(f, "%.15lf\n", x[i]);
		fclose(f);
	}

	//Ax = y умножение в разряженном формате
	//y - можно любой
	void MultMatrix(double *x, double *y, bool f)
	{
		if (f)
		{
			for (int i(0); i < n; i++)
				y[i] = diA[i] * x[i];


			for (int i(0); i < n; i++)
				for (int k = ig[i], m = ig[i + 1]; k < m; k++)
				{
					y[i] += gglA[k] * x[jg[k]];
					y[jg[k]] += gguA[k] * x[i];
				}
		}
		else
		{
			for (int i(0); i < n; i++)
				y[i] = diD[i] * x[i];
		}
	}

	void Solve()
	{
		MultMatrix(x, Ax, 1);
		for (int i(0); i < n; i++)
		{
			r[i] = pr[i] - Ax[i];
			z[i] = r[i];
		}
		MultMatrix(z, Ax, 1);
		swap(Ax, p);
		normar = scalar(r, r);
		for (iter = 0; iter<maxiter && abs(normar) > eps && abs(normar) > 1E-30; iter++)
		{
			double normap = scalar(p, p);
			double a = scalar(p, r) / normap;
			addmult(x, a, z, x);
			addmult(r, -a, p, r);
			MultMatrix(r, Ax, 1);
			double b = -scalar(p, Ax) / normap;
			addmult(r, b, z, z);
			normar -= a*a*normap;
			addmult(Ax, b, p, p);


		}
		int jgjtj = 0;
	}

	void FactorD()
	{
		try {
			for (int i(0); i < n; i++)
			{
				if (diA[i]<1E-30)
					throw - 1;
				diD[i] = 1. / sqrt(diA[i]);
			}
		}

		catch (int err)
		{
			printf_s("Divided by zero!");
		}
	}

	void SolveD()
	{
		MultMatrix(x, Ax, 1);
		for (int i(0); i < n; i++)
		{
			r[i] = diD[i] * (pr[i] - Ax[i]);
			z[i] = diD[i] * r[i];
		}
		MultMatrix(z, Ax, 1);
		MultMatrix(Ax, p, 0);
		normar = scalar(r, r);
		for (iter = 0; iter<maxiter && abs(normar) > eps && abs(normar) > 1E-30; iter++)
		{
			double normap = scalar(p, p);
			double a = scalar(p, r) / normap;
			addmult(x, a, z, x);
			addmult(r, -a, p, r);

			MultMatrix(r, Ax, 0);
			MultMatrix(Ax, Ay, 1);
			MultMatrix(Ay, Ay, 0);
			double b = -scalar(p, Ay) / normap;

			addmult(Ax, b, z, z);
			normar -= a*a*normap;
			addmult(Ay, b, p, p);

		}

	}

	void FactorLU()
	{
		int j0(0);
		for (int i(0); i < n; i++)
		{
			int m = ig[i + 1] - ig[i];
			double sum = 0;
			if (m != 0)
			{
				int jm = j0 + m;
				for (int j(j0); j < jm; j++)
				{
					int  stolbez(jg[j]);
					double sumL = 0, sumU = 0;
					if (j != j0)//npred!=0
						for (int l(ig[i]), nl(j), nu(ig[stolbez + 1]); l<nl; l++)
							for (int u(ig[stolbez]); u < nu; u++)
							{
								int gl = jg[l], gu = jg[u];
								if (gl == gu)
								{
									sumL += gguLU[u] * gglLU[l];
									sumU += gguLU[l] * gglLU[u];
								}
							}
					gglLU[j] = (gglA[j] - sumL) / diLU[jg[j]];
					gguLU[j] = (gguA[j] - sumU) / diLU[jg[j]];
					sum += gglLU[j] * gguLU[j];
				}
				j0 = jm;
			}
			diLU[i] = sqrt(diA[i] - sum);
		}
	}

	//Uy = x решение слау
	//н - любой
	void SLAEU(double *x, double *y)
	{
		for (int i(0); i < n; i++)y[i] = x[i];
		for (int i(n - 1); i >= 0; i--)
		{
			y[i] /= diLU[i];
			int kend = ig[i + 1];
			int k = ig[i];
			for (; k != kend; k++)
				y[jg[k]] -= y[i] * gguLU[k];
		}
	}

	//Ly = x решение слау
	//y - любой
	void SLAEL(double *x, double *y)
	{
		double sum;
		for (int i(0); i < n; i++)
		{
			sum = 0;
			int k0 = ig[i];
			int kend = ig[i + 1];
			for (int k = k0; k < kend; k++)
				sum += y[jg[k]] * gglLU[k];
			y[i] = (x[i] - sum) / diLU[i];
		}
	}

	void InitLU()
	{
		MultMatrix(x, r, 1);
		sub(r, pr, r);
		SLAEL(r, r);

		for (int i(0); i < n; i++)
			z[i] = r[i];

		SLAEU(z, z);

		MultMatrix(z, p, 1);
		SLAEL(p, p);

		normar = scalar(r, r);
	}

	void SolveLU()
	{
		InitLU();
		for (iter = 0; iter < maxiter && normar > eps && normar > 1E-30; iter++)
		{
			double normap = scalar(p, p);
			double a = scalar(p, r) / normap;
			addmult(x, a, z, x);
			addmult(r, -a, p, r);

			//(L-1)A(U-1)r = Ay
			//U-1r = Ax
			SLAEU(r, Ax);
			MultMatrix(Ax, Ay, 1);
			SLAEL(Ay, Ay);
			double b = -scalar(p, Ay) / normap;

			addmult(Ax, b, z, z);
			normar = abs(normar - a*a*normap);
			addmult(Ay, b, p, p);

			if (abs(a*a*normap) < 1E-20)
				InitLU();
		}
	}

	void ClearX()
	{
		for (int i(0); i < n; i++)
			x[i] = x0[i];
	}

	void Delete()
	{
		delete gglLU;
		delete gguLU;
		delete diLU;
		delete gglA;
		delete gguA;
		delete diA;
		delete diD;
		delete ig;
		delete jg;
		delete r;
		delete z;
		delete x;
		delete p;
		delete pr;
		delete Ax;
		delete Ay;
		delete x0;
	}
};