#include <iostream>
#include <vector>
#include "GMRESInterface.h"


using namespace std;
using namespace GMRESInterface;

int main()
{
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

    for (int i = 0; i < matrix_K.size(); i++)
    {
        for (int j = 0; j < matrix_K[i].size(); j++) {
            cout << matrix_K[i][j] << "\t";
        }
        cout << endl;
    }

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

    //GMRES
    for (int i = 0; i < 9; i++)
    {
        vec_q.push_back(0);
    }

    vector<vector<double>> temp_1 = { {1,2}, {1,2} };
    vector<vector<double>> temp_2 = { {3,4}, {3,4} };

    vector<vector<double>> temp_3 = MatrixByMatrix(temp_1, temp_2);
    vector<vector<double>> temp_4 = MatrixByDigit(temp_1, 1.0/2.0);
    for (int i = 0; i < temp_3.size(); i++)

    {
        for (int j = 0; j < temp_3[i].size(); j++) {
            cout << temp_3[i][j] << "\t";
        }
        cout << endl;
    }
    
    cout << endl << endl;
    for (int i = 0; i < temp_4.size(); i++)

    {
        for (int j = 0; j < temp_4[i].size(); j++) {
            cout << temp_4[i][j] << "\t";
        }
        cout << endl;
    }
    return 0;
}