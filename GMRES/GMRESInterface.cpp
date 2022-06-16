#include "GMRESInterface.h"

vector<vector<double>> GMRESInterface::Transponir(vector<vector<double>> matrix)
{
    size_t maxSize = 0;
    for (int i = 0; i < matrix.size(); i++)
    {
        maxSize = max(maxSize, matrix[i].size());
    }
    vector<vector<double>> temp(maxSize);
    // temp = NachPriblizh(temp);
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix[i].size(); j++)
        {
            temp[j].push_back(matrix[i][j]);
        }
    }
    return temp;
}

vector<vector<double>> GMRESInterface::NachPriblizh(vector<vector<double>> Matrix)
{
    vector<vector<double>> temp(Matrix.size());
    for (int i = 0; i < Matrix.size(); i++)
    {
        for (int j = 0; j < Matrix[i].size(); j++) temp[i].push_back(0.0);
    }
    return temp;
}

vector<double> GMRESInterface::NachPriblizh(vector<double> vec) {
    vector<double> temp(vec.size());
    for (int i = 0; i < vec.size(); i++)
    {
        temp[i] = 0.0;
    }
    return temp;
}

vector<vector<double>> GMRESInterface::RotateMatrix(vector<vector<double>> matrix_psi, vector<vector<double>> matrix_sigma, int ind)
{
    vector<vector<double>> temp(matrix_psi.size());
    for (int i = 0; i < matrix_psi.size(); i++)
    {
        for (int j = 0; j < matrix_psi.size(); j++) {
            if (i == j) temp[i].push_back(1.0);
            else temp[i].push_back(0.0);
        }
    }
    temp[ind][ind] = matrix_sigma[ind][ind] / (sqrtl(matrix_sigma[ind][ind] * matrix_sigma[ind][ind] + matrix_sigma[ind + 1][ind] * matrix_sigma[ind + 1][ind]));
    temp[ind + 1][ind + 1] = temp[ind][ind];

    temp[ind][ind + 1] = matrix_sigma[ind + 1][ind] / (sqrtl(matrix_sigma[ind][ind] * matrix_sigma[ind][ind] + matrix_sigma[ind + 1][ind] * matrix_sigma[ind + 1][ind]));
    temp[ind + 1][ind] = (-1.0) * temp[ind][ind + 1];
    //tarnspon
    /*temp[ind][ind] = matrix_sigma[ind][ind] / (sqrtl(matrix_sigma[ind][ind] * matrix_sigma[ind][ind] + matrix_sigma[ind][ind+1] * matrix_sigma[ind][ind+1]));
    temp[ind + 1][ind + 1] = temp[ind][ind];

    temp[ind+1][ind] = matrix_sigma[ind][ind+1] / (sqrtl(matrix_sigma[ind][ind] * matrix_sigma[ind][ind] + matrix_sigma[ind][ind+1] * matrix_sigma[ind][ind+1]));
    temp[ind][ind+1] = (-1.0) * temp[ind+1][ind];*/
    return temp;
}

vector<double> GMRESInterface::VecByMatrix(vector<double> Lvec, vector<vector<double>> RMatrix)
{
    vector<double> temp;
    for (int i = 0; i < Lvec.size(); i++)
    {
        double sum = 0.0;
        for (int j = 0; j < RMatrix.size(); j++)
        {
            sum += Lvec[j] * RMatrix[j][i];
        }
        temp.push_back(sum);
    }
    return temp;
    
}

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

vector<double> GMRESInterface::MatrixByVec(vector<vector<double>> LMatrix, vector<double> RVec)
{
    vector<double> temp;
    for (int i = 0; i < LMatrix.size(); i++)
    {
        double sum = 0.0;
        for (int j = 0; j < RVec.size(); j++)
        {
            sum += LMatrix[i][j] * RVec[j];
        }
        temp.push_back(sum);
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