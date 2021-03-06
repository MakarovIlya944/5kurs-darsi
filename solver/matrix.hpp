#pragma once
#include "net.hpp"
#include <set>
#include <map>

using namespace std;

struct Matrix
{
	int n;
	vector<int> ig, jg;
	vector<double> di, ggl, ggu, f;
};

class MatrixGenerator
{
	int *ig, *jg;
	double *di, *ggl, *ggu, *f;
	const double Beta[2]{1, 1};

	//Q[nt][nx]
	double **Q;

	//ѕересчитать! ст.227 соловейчик
	const double G[3][3] = {
		{2.33333333333333, -2.66666666666666, 0.33333333333333},
		{-2.66666666666666, 5.33333333333333, -2.66666666666666},
		{0.33333333333333, -2.66666666666666, 2.33333333333333}};
	const double GH[3][3] = {
		{0.50000000000000, -0.66666666666666, 0.16666666666666},
		{-0.66666666666666, 2.66666666666666, -2.00000000000000},
		{0.16666666666666, -2.00000000000000, 1.83333333333333}};
	const double M1[3][3] = {
		{0.1166666667, 0.06666666667, -0.01666666667},
		{0.06666666667, 0.2666666667, 0},
		{-0.01666666667, 0, 0.01666666667}};
	const double M2[3][3] = {
		{0.01666666667, 0., -0.01666666667},
		{0., 0.2666666667, 0.06666666667},
		{-0.01666666667, 0.06666666667, 0.1166666667}};
	const double MH1[3][3] = {
		{0.01190476190, 0.009523809524, -0.004761904762},
		{0.009523809524, 0.1142857143, 0.009523809524},
		{-0.004761904762, 0.009523809524, 0.01190476190}};
	const double MH2[3][3] = {
		{0.004761904762, -0.009523809524, -0.01190476190},
		{-0.009523809524, 0.1523809524, 0.05714285714},
		{-0.01190476190, 0.05714285714, 0.1047619048}};

	set<int> *edges;

	double LocalMatrixG[9][9]{}, LocalMatrixM[9][9]{}, LocalRightPart[9];
	double LocalMatrixC[9][9];
	double LocalMatrix[9][9];

	double LocalMatrixMr[3][3];
	double LocalMatrixMz[3][3];

	Net *net;
	map<double, Matrix> matricies;
	map<string, int> functions;

	int u(int i) { return i % 3; }

	int v(int i) { return i / 3; }

	void AccountNaturalCondition(int t);

	void AccountMainCondition(int t);

	void MakeLocalMatrixMassR(double r, double g1, double g2, double hr);

	void MakeLocalMatrixMassZ(double g1, double g2, double hz);

	void MakeLocalMatrix(int t, vector<int>::iterator elements);

public:
	MatrixGenerator() {}

	//Allocate memory for ig, init edges, ig - missing first 0; Q - answer
	//Generates ig, jg; allocate memor for ggl, ggu
	MatrixGenerator(Net *_net, double **&_Q);

	//elements - nvtr, GlobalPoints
	void CollectGlobalMatrix();

	void AccounConditions(int t);

	Matrix *GetMatrix(double t);

	map<double, Matrix> GetMatrix();

	void SaveMatrix();

	void ClearMatrix();

	~MatrixGenerator();
};
