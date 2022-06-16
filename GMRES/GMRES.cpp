#include <iostream>
#include <fstream>
#include <vector>
#include "GMRESInterface.h"
#include "Data.h"

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
    //int N = 4;
    int m = 5;
    double eps = 10E-3;
    vector<double> vec_q;
    vector<vector<double>> matrix_K((data::N - 1) * (data::N - 1));
    //vectorPrintFile(matrix_K, fout);
    vector<double> vec_X((data::N - 1) * (data::N - 1));
    

    matrix_K = data::fillingVectorK(data::N);

    //MatrixView(matrix_K);

    cout << endl << endl;

    vec_X = data::fillVectorX(data::N);

    //MatrixView(vec_X);

    vector<vector<double>> EMatrix((data::N - 1) * (data::N - 1));

    for (int i = 0; i < (data::N - 1) * (data::N - 1); i++)
    {
        for (int j = 0; j < (data::N - 1) * (data::N - 1); j++)
        {
            if (i == j) EMatrix[i].push_back(1);
            else EMatrix[i].push_back(0);
        }
    }


    vec_q = vec_X;
    vector<double> vec_R_0 = vec_X;
    vector<double> vec_R_1;
    vec_R_0 = VecMinusVec(vec_X, MatrixByVec(matrix_K, vec_q));
    vector<double> vec_R_ = vec_X;
    //1 stage

    vec_q = NachPriblizh(vec_q);
    // EuqlidNorm(MatrixByVec(EMatrix, vec_R_0)) <= eps * EuqlidNorm(MatrixByVec(EMatrix, vec_R_0))
    do {
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
        sigma = Transponir(sigma);

        //vectorPrintFile(teta, fout);
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
            sigma = MatrixByMatrix(matrix_psi, sigma);
            //vectorPrintFile(sigma, fout);
            vec_e_1_sh = MatrixByVec(matrix_psi, vec_e_1_sh);
            //vectorPrintFile(vec_e_1_sh, fout);
        }
        sigma.pop_back();
        //vectorPrintFile(sigma, fout);
        //vectorPrintFile(vec_e_1_sh, fout);
        vector<double> vec_z(m+1);
        vec_z = NachPriblizh(vec_z);

        for (int i = 0; i < m; i++)
        {
            double summm = 0.0;
            for (int j = 0; j < i; j++)
            {
                summm += sigma[m-1-i][m-1-j]*vec_z[m-1-j];
            }
            vec_z[m - i] = (vec_e_1_sh[m - i] - summm) / sigma[m-1 - i][m-1 - i];
        }
        //vectorPrintFile(vec_z, fout);
        vector<double> temp(teta[0].size());
        for (int i = 0; i < m; i++)
        {
            temp = VecAddVec(temp, VecByDigit(teta[i], vec_z[i]));
        }
        vec_q = VecAddVec(vec_q, temp);
        vectorPrintFile(vec_q, fout);
        
    } while (EuqlidNorm(MatrixByVec(EMatrix, vec_R_)) <= eps * EuqlidNorm(MatrixByVec(EMatrix, vec_R_0)));

    //vectorPrintFile(MatrixByVec(matrix_K, vec_q), fout);
    fout.close();
    return 0;
}