#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "GMRESInterface.h"
#include "Data.h"

using namespace std;
using namespace GMRESInterface;



void MatrixView(vector<vector<double>> &matrix) {
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix[i].size(); j++) {
            cout << matrix[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl << endl;
}

void MatrixView(vector<double> &vec) {
    for (int j = 0; j < vec.size(); j++) {
        cout << vec[j] << endl;
    }
    cout << endl << endl;
}

void vectorPrintFile(vector<vector<double>> &vec, ostream& fout) {
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

void vectorPrintFile(vector<double> &vec, ostream& fout) {
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
    int m = 5;
    double eps = 1E-3;

    string fileName = "N";
    fileName.append(std::to_string(data::N));
    fileName += "_eps";
    fileName.append(std::to_string(eps));
    fileName += "_m";
    fileName.append(std::to_string(m));
    fileName += ".txt";
    ofstream fout(fileName);

    vector<double> vec_q((data::N - 1) * (data::N - 1), 0.0);
    vector<vector<double>> matrix_K((data::N - 1) * (data::N - 1));
    
    vector<double> vec_X((data::N - 1) * (data::N - 1));
    

    matrix_K = data::fillingVectorK(matrix_K, data::N);
    //vectorPrintFile(matrix_K, fout);
    vec_X = data::fillVectorX(vec_X, data::N);
    //vectorPrintFile(vec_X, fout);

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
    //vector<vector<double>> A = { {4, 3, 1}, {3, 5 ,7}, {4, 2, 6} };
    vector<vector<double>> L = matrix_K;
    vector<vector<double>> U = matrix_K;

    L = NachPriblizh(L);
    U = NachPriblizh(U);

    RazlozhenieLU(matrix_K, L, U);

    //vector<double> rv = { 3, 2, 1 };
    
    //vector<double> temp_ = Pereobuslav(L, U, rv);


    //vec_q = vec_X;
    vec_q = NachPriblizh(vec_q);
    vector<double> vec_R_0 = vec_X;
    vector<double> vec_R_1;
    vector<double> TEMP = MatrixByVec(matrix_K, vec_q);
    vec_R_0 = VecMinusVec(vec_X, TEMP);
    vector<double> vec_R_ = vec_R_0;
    //1 stage
    TEMP = MatrixByVec(EMatrix, vec_R_);
    long double NormTemp = EuqlidNorm(TEMP);
    TEMP = MatrixByVec(EMatrix, vec_R_0);
    long double NormZero = eps * EuqlidNorm(TEMP);
    /*cout << EuqlidNorm(MatrixByVec(EMatrix, vec_R_)) << endl;
    cout << eps * EuqlidNorm(MatrixByVec(EMatrix, vec_R_0));*/
    // EuqlidNorm(MatrixByVec(EMatrix, vec_R_0)) <= eps * EuqlidNorm(MatrixByVec(EMatrix, vec_R_0))
    int counter = 0;
    while (NormTemp > NormZero) {
        //2 stage
        counter++;
        TEMP = MatrixByVec(matrix_K, vec_q);
        vec_R_ = VecMinusVec(vec_X, TEMP);
        NormTemp = EuqlidNorm(vec_R_);
        vec_R_1 = MatrixByVec(EMatrix, vec_R_);
        //3 stage
        long double varrho = EuqlidNorm(vec_R_1);
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
            //vectorPrintFile(rho, fout);
            //vector<vector<double>> sigma(m);
            vector<double> _sigma(m + 1);
            _sigma = NachPriblizh(_sigma);
            int i;
            for (i = 0; i <= j; i++)
            {
                //8 stage
                _sigma[i] = ScalarVecByVec(rho, teta[i]);
                //9 stage
                TEMP = VecByDigit(teta[i], _sigma[i]);
                rho = VecMinusVec(rho, TEMP);
                //vectorPrintFile(rho, fout);
            } //10 stage
            //11 stage
            _sigma[i] = EuqlidNorm(rho);
            sigma.push_back(_sigma);
            if (EuqlidNorm(rho) == 0) break;
            //12 stage
            else teta.push_back(VecByDigit(rho, 1.0 / EuqlidNorm(rho)));
        } //13 stage
        sigma = Transponir(sigma);
        //vectorPrintFile(sigma, fout);
        //vectorPrintFile(teta, fout);
        vector<double> vec_e_1(m+1);
        vec_e_1 = NachPriblizh(vec_e_1);
        vec_e_1.front() = 1.0;
        vector<double> vec_e_1_sh(m+1);
        vector<vector<double>> matrix_psi(m + 1);
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
        vec_e_1_sh.pop_back();
        //vectorPrintFile(sigma, fout);
        //vectorPrintFile(vec_e_1_sh, fout);
        vector<double> vec_z(m);
        vec_z = NachPriblizh(vec_z);

        for (int i = 0; i < m; i++)
        {
            double summm = 0.0;
            for (int j = 0; j < i; j++)
            {
                summm += sigma[m - 1 - i][m - 1 - j] * vec_z[m - 1 - j];
            }
            vec_z[m - 1 - i] = (vec_e_1_sh[m - 1 - i] - summm) / sigma[m - 1 - i][m - 1 - i];
        }
        //vectorPrintFile(vec_z, fout);
        vector<double> temp(teta[0].size(), 0.0);
        for (int i = 0; i < m; i++)
        {
            TEMP = VecByDigit(teta[i], vec_z[i]);
            temp = VecAddVec(temp, TEMP);
        }
        vec_q = VecAddVec(vec_q, temp);
        //vectorPrintFile(vec_q, fout);
        cout << "Iteration #" << counter - 1 << " end!\n\n";
        fout << "Iteration #" << counter - 1 << " end!\t" << NormTemp << endl << endl;

    };

    vector<double>temp_uzl;

    for (int i = 1; i < data::N; i++)
    {
        for (int j = 1; j < data::N; j++) {
            double x = data::x(j, data::h);
            double y = data::y(i, data::h);
            temp_uzl.push_back(x * x * y * y - y * y * x - x * x * y + x * y);
            //temp_uzl[i][j] = ((double)i * (double)i * h * h - (double)i * h) * ((double)j * (double)j * h * h - (double)j * h);
        }
    }

    //vectorPrintFile(vec_q, fout);
    //vectorPrintFile(vec_X, fout);
    //vectorPrintFile(temp_uzl, fout);

    if (vec_q.size() == temp_uzl.size()) cout << "size equal";

    else cout << "size not equal";
    double max = 0.0;
    int index_max = 0;
    for (int i = 0; i < (data::N-1)*(data::N-1); i++)
    {
        if (max < abs(vec_q[i] - temp_uzl[i])) {
            max = abs(vec_q[i] - temp_uzl[i]);
            index_max = i;
        }
    }

    fout << endl << endl << max << "\t" << index_max;

    
    fout.close();
    return 0;
}