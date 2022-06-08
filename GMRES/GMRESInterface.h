#pragma once
#include <vector>
#include <math.h>
using namespace std;

namespace GMRESInterface
{
	vector<vector<double>> MatrixByMatrix(vector<vector<double>> LMatrix, vector<vector<double>> RMatrix);
	vector<vector<double>> MatrixByDigit(vector<vector<double>> LMatrix, double digit);
	vector<double> VecByDigit(vector<double> vec, double digit);
	vector<double> VecMinusVec(vector<double> LVec, vector<double> RVec);
	vector<double> VecAddVec(vector<double> LVec, vector<double> RVec);
	long double ScalarVecByVec(vector<double> LVec, vector<double> RVec);
	long double EuqlidNorm(vector<vector<double>> Matrix);
	long double EuqlidNorm(vector<double> vec);

};

