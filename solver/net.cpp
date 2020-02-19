#include "net.hpp"
#include <fstream>
#include <iostream>
#include <list>
#include <iterator>

using namespace std;

Net* NetGenerator::Generate(vector<double> R, vector<double> Z, double step, vector<double> r_steps, 
	vector<double> z_steps, vector<int> W, vector<Pointi> area, vector<int> B, vector<int> V, vector<Pointi> border,
	vector<double> T, double t_step, vector<double> t_steps)
{
	// materials.size() // 2 == 0
	// r_steps.size() // 2 == 0 && r_steps.size() / 2 == R.size() - 1
	// Z.size() // 2 == 0 && z_steps.size() / 2 == Z.size() - 1

	Net *net = new Net();
	try
	{
		int nr = R.size(), nz = Z.size();
		step *= 0.5;
		//double* kr = new double[nr - 1], * kz = new double[nz - 1];
		//int* hr = new int[nr - 1], nhr(0), * hz = new int[nz - 1], nhz(0);
		//r_steps ~ hr[0], kr[0], hr[1], kr[1], ...
		int nhr(0), nhz(0);
		for (auto i = r_steps.begin(); i != r_steps.end(); i+=2)
		{
			*i *= 2;
			nhr += (int)(*i);
		}
		for (auto i = z_steps.begin(); i != z_steps.end(); i += 2)
		{
			*i *= 2;
			nhz += (int)(*i);
		}

		int nw = W.size(), i(0), j(0), k(0);
		double* Rh = new double[nhr + 1], * Zh = new double[nhz + 1];
		for (int i(0), k(0); i < nr - 1; i++)
		{
			for (int j(0); j < r_steps[i*2]; j++, k++)
				Rh[k] = R[i] + j * step * r_steps[i*2+1];
			r_steps[i*2] += i != 0 ? r_steps[(i - 1)*2] : 0;
		}
		Rh[nhr] = R[nr - 1];

		for (int i(0), k(0); i < nz - 1; i++)
		{
			for (int j(0); j < z_steps[i*2]; j++, k++)
				Zh[k] = Z[i] + j * step * z_steps[i*2+1];
			z_steps[i*2] += i != 0 ? z_steps[(i - 1)*2] : 0;
		}
		Zh[nhz] = Z[nz - 1];

		net->rNum = nhr + 1;
		net->zNum = nhz + 1;
		net->Num = net->rNum * net->zNum;
		net->GlobalNet = vector<Pointd>(net->Num) ;// new Pointd[net->Num];
		int n = nhr * nhz / 4;
		net->nvtr = vector<int>(n * 9);// new int[n * 9];
		net->nvcat = vector<int>(n);//new int[n];
		net->NumEl = n;
		double tmp;
		for (int i(0), w(0), currentelem(0); i < net->zNum; i++)
		{
			tmp = Zh[i];
			for (int j(0); j < net->rNum; j++)
			{
				if ((i % 2) && (j % 2))
				{
					for (int k(0), n(nw * 2); k < n && w == 0; k += 2)
					{
						if ((area[k].r == 0 || r_steps[(area[k].r - 1)*2] < j) && j < r_steps[(area[k + 1].r - 1)*2])
							if (i < z_steps[(area[k + 1].z - 1)*2] && (area[k].z == 0 || z_steps[(area[k].z - 1)*2] < i))
								w = W[k / 2];
					}
					net->nvcat[currentelem++] = w;
					w = 0;
				}
				net->GlobalNet[i * net->rNum + j] = Pointd(Rh[j], tmp);
			}
		}

		GenerationNVTR(net);

		net->nvr1.clear();
		net->nvr2.clear();
		net->nvr3.clear();

		int nb = B.size();

		for (int i(0), j(0); i < nb; i++, j++)
		{
			list<int>* iterator = &(net->nvr1);

			switch (B[i])
			{
			case 2:
				iterator = &(net->nvr2);
				break;
			case 3:
				iterator = &(net->nvr3);
			}
			iterator->push_back(V[i]);
			if (border[j].r == border[j + 1].r)
			{
				int l = border[j].z != 0 ? z_steps[(border[j].z - 1)*2] : 0;
				int k = border[j].r != 0 ? r_steps[(border[j].r - 1)*2] : 0;
				if (border[j].z == border[j + 1].z)
					iterator->push_back(l * net->rNum + k);
				int nl(z_steps[(border[++j].z - 1)*2]);
				for (; l <= nl; l++)
					iterator->push_back(l * net->rNum + k);
			}
			else
			{
				int l = border[j].r != 0 ? r_steps[(border[j].r - 1)*2] : 0;
				int k = border[j].z != 0 ? z_steps[(border[j].z - 1)*2] : 0;
				int nl(r_steps[(border[++j].r - 1)*2]);
				for (; l <= nl; l++)
					iterator->push_back(k * net->rNum + l);
			}
			iterator->push_back(-1);
		}

		int nt = T.size(), nht(0);
		for (auto i= t_steps.begin(); i != t_steps.end(); i+=2)
			nht += (int)(*i);

		net->GlobalTime = vector<double>(nht + 1);// new double[nht + 1];
		net->GlobalTime[0] = T[0];
		for (int i(0), k(0); k < nht; i++)
		{
			tmp = t_step * t_steps[i*2+1];
			for (int j(0); j < t_steps[i*2]; j++, k++)
				net->GlobalTime[k + 1] = net->GlobalTime[k] + tmp;
		}
		net->NumT = nht + 1;

		delete[] Rh;
		delete[] Zh;
	}
	catch (const std::exception & exception)
	{
		printf_s("Error: %s", exception.what());
	}
	if (isSavingNet)
	{
		ofstream out;
		out.open(fileName);

		for (int i(0); i < net->NumT; i++)
			out << net->GlobalTime[i];
		for (int i(0); i < net->rNum; i++)
			out << net->GlobalNet[i].r;
		for (int i(0), j(0); i < net->zNum; i++, j += net->rNum)
			out << net->GlobalNet[i].z;

		out.close();
	}

	return net;
}


Net* NetGenerator::GenerateFromFiles(const char* FNameN, const char* FNameB, const char* FNameT)
{
	Net *net = new Net();
	try
	{
		ifstream in(FNameN);
		int nr, nz;
		in >> nr;
		if (isLoging)
			printf_s("nr: %d\n", nr);
		double* R = new double[nr];
		for (int i(0); i < nr; i++)
		{
			in >> R[i];
			if (isLoging)
				printf_s("R[%d]: %f\n", i, R[i]);
		}

		in >> nz;
		if (isLoging)
			printf_s("nz: %d\n", nz);
		double* Z = new double[nz];
		for (int i(0); i < nz; i++)
		{
			in >> Z[i];
			if (isLoging)
				printf_s("Z[%d]: %f\n", i, Z[i]);
		}

		double step;
		in >> step;
		if (isLoging)
			printf_s("step: %f\n", step);
		step *= 0.5;
		double* kr = new double[nr - 1], * kz = new double[nz - 1];
		int* hr = new int[nr - 1], nhr(0), * hz = new int[nz - 1], nhz(0);
		for (int i(0); i < nr - 1; i++)
		{
			in >> hr[i];
			in >> kr[i];
			if (isLoging)
				printf_s("hr[%d]: %d\nkr[%d]: %f\n", i, hr[i], i,  kr[i]);
			hr[i] *= 2;
			nhr += hr[i];
		}
		for (int i(0); i < nz - 1; i++)
		{
			in >> hz[i];
			in >> kz[i];
			if (isLoging)
				printf_s("hz[%d]: %d\nkz[%d]: %f\n", i, hz[i], i, kz[i]);
			hz[i] *= 2;
			nhz += hz[i];
		}

		int nw;
		in >> nw;
		if (isLoging)
			printf_s("step: %d\n", nw);
		int* W = new int[nw];
		Pointi* area = new Pointi[nw * 2]{};
		for (int i(0), j(0); i < nw; i++)
		{
			in >> W[i];
			in >> area[j].r;
			in >> area[j++].z;
			in >> area[j].r;
			in >> area[j++].z;
			if (isLoging)
				printf_s("W[%d]: %d\nArea[%d]: [%d %d - %d %d]", i, W[i], (j-1) / 2, area[j-2].r, area[j-2].z, area[j-1].r, area[j-1].z);
		}

		double* Rh = new double[nhr + 1], * Zh = new double[nhz + 1];
		for (int i(0), k(0); i < nr - 1; i++)
		{
			for (int j(0); j < hr[i]; j++, k++)
				Rh[k] = R[i] + j * step * kr[i];
			hr[i] += i != 0 ? hr[i - 1] : 0;
		}
		Rh[nhr] = R[nr - 1];

		for (int i(0), k(0); i < nz - 1; i++)
		{
			for (int j(0); j < hz[i]; j++, k++)
				Zh[k] = Z[i] + j * step * kz[i];
			hz[i] += i != 0 ? hz[i - 1] : 0;
		}
		Zh[nhz] = Z[nz - 1];

		delete[] R;
		delete[] Z;
		net->rNum = nhr + 1;
		net->zNum = nhz + 1;
		net->Num = net->rNum * net->zNum;
		net->GlobalNet = vector<Pointd>(net->Num);// = new Pointd[net->Num];
		int n = nhr * nhz / 4;
		net->nvtr = vector<int>(n * 9);// = new int[n * 9];
		net->nvcat = vector<int>(n);//= new int[n];
		net->NumEl = n;
		double tmp;
		for (int i(0), w(0), currentelem(0); i < net->zNum; i++)
		{
			tmp = Zh[i];
			for (int j(0); j < net->rNum; j++)
			{
				if ((i % 2) && (j % 2))
				{
					for (int k(0), n(nw * 2); k < n && w == 0; k += 2)
					{
						if ((area[k].r == 0 || hr[area[k].r - 1] < j) && j < hr[area[k + 1].r - 1])
							if (i < hz[area[k + 1].z - 1] && (area[k].z == 0 || hz[area[k].z - 1] < i))
								w = W[k / 2];
					}
					net->nvcat[currentelem++] = w;
					w = 0;
				}
				net->GlobalNet[i * net->rNum + j] = Pointd(Rh[j], tmp);
			}
		}

		GenerationNVTR(net);
		in.close();

		in = ifstream(FNameB);

		net->nvr1.clear();
		net->nvr2.clear();
		net->nvr3.clear();

		int nb;
		in >> nb;
		int* B = new int[nb], * V = new int[nb];
		Pointi* border = new Pointi[nb * 2]{};
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
			list<int>* iterator = &(net->nvr1);

			switch (B[i])
			{
			case 2:
				iterator = &(net->nvr2);
				break;
			case 3:
				iterator = &(net->nvr3);
			}
			iterator->push_back(V[i]);
			if (border[j].r == border[j + 1].r)
			{
				int l = border[j].z != 0 ? hz[border[j].z - 1] : 0;
				int k = border[j].r != 0 ? hr[border[j].r - 1] : 0;
				if (border[j].z == border[j + 1].z)
					iterator->push_back(l * net->rNum + k);
				int nl(hz[border[++j].z - 1]);
				for (; l <= nl; l++)
					iterator->push_back(l * net->rNum + k);
			}
			else
			{
				int l = border[j].r != 0 ? hr[border[j].r - 1] : 0;
				int k = border[j].z != 0 ? hz[border[j].z - 1] : 0;
				int nl(hr[border[++j].r - 1]);
				for (; l <= nl; l++)
					iterator->push_back(k * net->rNum + l);
			}
			iterator->push_back(-1);
		}

		in.close();
		in.open(FNameT);

		int nt;
		in >> nt;
		double* T = new double[nt];
		for (int i(0); i < nt; i++)
			in >> T[i];

		in >> step;

		double* kt = new double[nt - 1];
		int* ht = new int[nt - 1], nht(0);
		for (int i(0); i < nt - 1; i++)
		{
			in >> ht[i];
			in >> kt[i];
			nht += ht[i];
		}

		net->GlobalTime = vector<double>(nht + 1);// = new double[nht + 1];
		net->GlobalTime[0] = T[0];
		for (int i(0), k(0); k < nht; i++)
		{
			tmp = step * kt[i];
			for (int j(0); j < ht[i]; j++, k++)
				net->GlobalTime[k + 1] = net->GlobalTime[k] + tmp;
		}
		net->NumT = nht + 1;
		in.close();


		delete[] kt;
		delete[] ht;
		delete[] T;
		delete[] Rh;
		delete[] Zh;
		delete[] hz;
		delete[] hr;
		delete[] kr;
		delete[] kz;
		delete[] W;
		delete[] border;

	}
	catch (const std::exception& exception)
	{
		printf_s("Error: %s", exception.what());
	}
	if (isSavingNet)
	{
		ofstream out;
		out.open(fileName);

		for (int i(0); i < net->NumT; i++)
			out << net->GlobalTime[i];
		for (int i(0); i < net->rNum; i++)
			out << net->GlobalNet[i].r;
		for (int i(0), j(0); i < net->zNum; i++, j += net->rNum)
			out << net->GlobalNet[i].z;

		out.close();
	}

	return net;
}

void NetGenerator::GenerationNVTR(Net *net)
{
	int zn = net->zNum - 1, rn = net->rNum - 1;
	for (int z(0), current(0); z < zn; z += 2)
		for (int r(0); r < rn; r += 2)
		{
			int index = r + z * net->rNum;
			net->nvtr[current++] = index;
			net->nvtr[current++] = index + 1;
			net->nvtr[current++] = index + 2;

			net->nvtr[current++] = index + net->rNum;
			net->nvtr[current++] = index + net->rNum + 1;
			net->nvtr[current++] = index + net->rNum + 2;

			net->nvtr[current++] = index + 2 * net->rNum;
			net->nvtr[current++] = index + 2 * net->rNum + 1;
			net->nvtr[current++] = index + 2 * net->rNum + 2;
		}
}

	Net::~Net()
	{
		GlobalNet.clear();
		nvcat.clear();
		nvtr.clear();
		nvr1.clear();
		nvr2.clear();
		nvr3.clear();
	}
