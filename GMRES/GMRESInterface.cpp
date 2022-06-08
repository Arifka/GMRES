#include "GMRESInterface.h"

vector<vector<double>> GMRESInterface::MatrixByMatrix(vector<vector<double>> LMatrix, vector<vector<double>> RMatrix)
{
    vector<vector<double>> temp(LMatrix.size());
    for (int i = 0; i < LMatrix.size(); i++)
    {
        for (int j = 0; j < RMatrix[i].size(); j++)
        {
            double sum = 0.0;
            for (int k = 0; k < RMatrix.size(); k++)
            {
                sum += LMatrix[i][k] * RMatrix[k][j];
            }
            temp[i].push_back(sum);
        }
    }
    return temp;
}

vector<vector<double>> GMRESInterface::MatrixByDigit(vector<vector<double>> LMatrix, double digit)
{
    vector<vector<double>> temp = LMatrix;
    for (int i = 0; i < LMatrix.size(); i++)
    {
        for (int j = 0; j < LMatrix[i].size(); j++) {
            temp[i][j] *= digit;
        }
    }
    return temp;
}

vector<double> GMRESInterface::VecByDigit(vector<double> vec, double digit)
{
    vector<double> temp = vec;
    for (int i = 0; i < vec.size(); i++)
    {
        temp[i] *= digit;
    }
    return temp;
}

vector<double> GMRESInterface::VecMinusVec(vector<double> LVec, vector<double> RVec)
{
    vector<double> temp(LVec.size());
    for (int i = 0; i < temp.size(); i++)
    {
        temp[i] = LVec[i] - RVec[i];
    }
    return temp;
}

vector<double> GMRESInterface::VecAddVec(vector<double> LVec, vector<double> RVec)
{
    vector<double> temp(LVec.size());
    for (int i = 0; i < temp.size(); i++)
    {
        temp[i] = LVec[i] + RVec[i];
    }
    return temp;
}

long double GMRESInterface::ScalarVecByVec(vector<double> LVec, vector<double> RVec)
{
    long double sum = 0.0;
    for (int i = 0; i < LVec.size(); i++)
    {
        sum += LVec[i] * RVec[i];
    }
    return sum;
}

long double GMRESInterface::EuqlidNorm(vector<vector<double>> Matrix)
{
    long double sum = 0.0;
    for (int i = 0; i < Matrix.size(); i++)
    {
        for (int j = 0; j < Matrix[i].size(); j++)
        {
            sum += Matrix[i][j] * Matrix[i][j];
        }
    }
    return sqrtl(sum);
}

long double GMRESInterface::EuqlidNorm(vector<double> vec) 
{
    long double sum = 0.0;
    for (int i = 0; i < vec.size(); i++)
    {
        sum += vec[i] * vec[i];
    }
    return sqrtl(sum);
}