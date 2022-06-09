#pragma once
#include <vector>
#include <math.h>
using namespace std;

namespace GMRESInterface
{
	vector<vector<double>> NachPriblizh(vector<vector<double>> Matrix);
	vector<double> NachPriblizh(vector<double> vec);
	vector<vector<double>> MatrixByMatrix(vector<vector<double>> LMatrix, vector<vector<double>> RMatrix);
	vector<double> MatrixByVec(vector<vector<double>> LMatrix, vector<double> RVec);
	vector<vector<double>> MatrixByDigit(vector<vector<double>> LMatrix, double digit);
	vector<double> VecByDigit(vector<double> vec, double digit);
	vector<double> VecMinusVec(vector<double> LVec, vector<double> RVec);
	vector<double> VecAddVec(vector<double> LVec, vector<double> RVec);
	long double ScalarVecByVec(vector<double> LVec, vector<double> RVec);
	long double EuqlidNorm(vector<vector<double>> Matrix);
	long double EuqlidNorm(vector<double> vec);

};

