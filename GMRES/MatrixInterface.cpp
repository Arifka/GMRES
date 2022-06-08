#include "MatrixInterface.h"

vector<vector<double>> MatrixInterface::MatrixMultiply(vector<vector<double>> LMatrix, vector<vector<double>> RMatrix)
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
