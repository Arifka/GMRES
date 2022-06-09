#include <iostream>
#include <vector>
#include "GMRESInterface.h"

using namespace std;
using namespace GMRESInterface;

void MatrixView(vector<vector<double>> matrix) {
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix[i].size(); j++) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl << endl;
}

void MatrixView(vector<double> vec) {
    for (int j = 0; j < vec.size(); j++) {
        cout << vec[j] << endl;
    }
    cout << endl << endl;
}

int main()
{
    int m = 10;
    double eps = 10E-3;
    vector<double> vec_q;
    vector<vector<double>> matrix_K(9);
    vector<double> vec_X;
    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++) {
            if (i == j) {
                matrix_K[i].push_back(4.0);
            }
            else if ((abs(j - i) == 1 or abs(j - i) == 3)) matrix_K[i].push_back(-1.0);
            else
            {
                matrix_K[i].push_back(0);
            }
        }
    }

    matrix_K[2][3] = 0;
    matrix_K[3][2] = 0;
    matrix_K[5][6] = 0;
    matrix_K[6][5] = 0;

    MatrixView(matrix_K);

    cout << endl << endl;

    vec_X.push_back(1.0 / 24.0);
    vec_X.push_back(41.0 / 768.0);
    vec_X.push_back(1.0 / 24.0);
    vec_X.push_back(41.0 / 768.0);
    vec_X.push_back(13.0 / 192.0);
    vec_X.push_back(41.0 / 768.0);
    vec_X.push_back(1.0 / 24.0);
    vec_X.push_back(41.0 / 768.0);
    vec_X.push_back(1.0 / 24.0);

    vector<vector<double>> EMatrix(9);

    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            if (i == j) EMatrix[i].push_back(1);
            else EMatrix[i].push_back(0);
        }
    }

    vector<vector<double>> temp_1 = { { 1, 2}, { 5, 3} };
    vector<double> temp_2 = { 7, 12 };

    vector<double> temp_3 = MatrixByVec(temp_1, temp_2);

    MatrixView(MatrixByVec(temp_1, temp_2));

    vector<double> vec_R_0 = vec_X;
    vector<double> vec_R_1;
    
    while (EuqlidNorm(MatrixByVec(EMatrix, vec_R_0)) <= eps * EuqlidNorm(MatrixByVec(EMatrix, vec_R_0)))
    {
        //1 stage
        vec_q = NachPriblizh(vec_q);
        //2 stage
        vec_R_0 = VecMinusVec(vec_X, MatrixByVec(matrix_K, vec_q));
        vec_R_1 = MatrixByVec(EMatrix, vec_R_0);
        //3 stage
        double varrho = EuqlidNorm(vec_R_1);
        //4 stage
        vector<vector<double>> teta;
        teta.push_back(VecByDigit(vec_R_1, 1.0/varrho));
        //5 stage

    }

    //GMRES
    for (int i = 0; i < 9; i++)
    {
        vec_q.push_back(0);
    }
    
    return 0;
}