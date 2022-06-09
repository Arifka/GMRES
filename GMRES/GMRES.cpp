#include <iostream>
#include <fstream>
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

void vectorPrintFile(vector<vector<double>> vec, ostream& fout) {
    fout.setf(ios::left);
    for (int i = 0; i < vec.size(); i++) {
        for (int j = 0; j < vec[i].size(); j++) {
            fout.width(20);
            fout << vec[i][j];
        }
        fout << endl;
    }
    fout.unsetf(ios::left);
    fout << endl << endl;
}

void vectorPrintFile(vector<double> vec, ostream& fout) {
    fout.setf(ios::left); 
        for (int j = 0; j < vec.size(); j++) {
            fout << vec[j];
            fout << endl;
        }
    fout.unsetf(ios::left);
    fout << endl << endl;
}

int main()
{
    ofstream fout("output.txt");
    int m = 8;
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

    vec_q = vec_X;
    vector<double> vec_R_0 = vec_X;
    vector<double> vec_R_1;
    vec_R_0 = VecMinusVec(vec_X, MatrixByVec(matrix_K, vec_q));
    vector<double> vec_R_ = vec_X;
    // EuqlidNorm(MatrixByVec(EMatrix, vec_R_0)) <= eps * EuqlidNorm(MatrixByVec(EMatrix, vec_R_0))
    do {
        //1 stage
        vec_q = NachPriblizh(vec_q);
        //2 stage
        vec_R_ = VecMinusVec(vec_X, MatrixByVec(matrix_K, vec_q));
        vec_R_1 = MatrixByVec(EMatrix, vec_R_);
        //3 stage
        double varrho = EuqlidNorm(vec_R_1);
        //4 stage
        vector<vector<double>> teta;
        teta.push_back(VecByDigit(vec_R_1, 1.0 / varrho));
        //5 stage
        vector<vector<double>> sigma;
        for (int j = 0; j < m; j++)
        {
            //6 stage
            vector<double> rho;
            vector<double> rho_ = MatrixByVec(matrix_K, teta[j]);
            rho = MatrixByVec(EMatrix, rho_);
            //7 stage
            //vector<vector<double>> sigma(m);
            vector<double> _sigma(m+1);
            _sigma = NachPriblizh(_sigma);
            int i;
            for (i = 0; i <= j; i++)
            {
                //8 stage
                _sigma[i] = ScalarVecByVec(rho, teta[i]);
                //9 stage
                rho = VecMinusVec(rho, VecByDigit(teta[i], _sigma[i]));
            } //10 stage
            //11 stage
            _sigma[i] = EuqlidNorm(rho);
            sigma.push_back(_sigma);
            if (EuqlidNorm(rho) == 0) break;
            //12 stage
            else teta.push_back(VecByDigit(rho, 1.0 / EuqlidNorm(rho)));
        } //13 stage
        vectorPrintFile(sigma, fout);
        vector<double> vec_e_1(m+1);
        vec_e_1 = NachPriblizh(vec_e_1);
        vec_e_1.front() = 1.0;
        vector<double> vec_e_1_sh(m+1);
        vector<vector<double>> matrix_psi(m+1);
        vec_e_1_sh = VecByDigit(vec_e_1, varrho);
        for (int i = 0; i < m; i++)
        {
            matrix_psi = RotateMatrix(matrix_psi, sigma, i);
            //vectorPrintFile(matrix_psi, fout);
            sigma = MatrixByMatrix(sigma, matrix_psi);
            //vectorPrintFile(sigma, fout);
            vec_e_1_sh = VecByMatrix(vec_e_1_sh, matrix_psi);
            //vectorPrintFile(vec_e_1_sh, fout);
        }
        for (int i = 0; i < m; i++)
        {
            sigma[i].pop_back();
        }
        vectorPrintFile(sigma, fout);
        vectorPrintFile(vec_e_1_sh, fout);
        vector<double> vec_z(m+1);
        vec_z = NachPriblizh(vec_z);

        for (int i = 0; i < m; i++)
        {
            double summm = 0.0;
            for (int j = 0; j < i; j++)
            {
                summm += sigma[j][i]*vec_z[m-j];
            }
            vec_z[m - i] = (vec_e_1_sh[m - i] - summm) / sigma[m-1 - i][m-1 - i];
        }
        vectorPrintFile(vec_z, fout);
        vector<double> temp(teta[0].size());
        for (int i = 0; i < m; i++)
        {
            temp = VecAddVec(temp, VecByDigit(teta[i], vec_z[i]));
        }
        vec_q = VecAddVec(vec_q, temp);
        vectorPrintFile(vec_q, fout);
        
    } while (EuqlidNorm(MatrixByVec(EMatrix, vec_R_)) <= eps * EuqlidNorm(MatrixByVec(EMatrix, vec_R_0)));

    fout.close();
    return 0;
}