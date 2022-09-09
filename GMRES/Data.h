#pragma once
#include <vector>


namespace data {
	const int N = 4;
	const double l = 1.0;
	const double h = l / static_cast<double>(N);
    
    double x(int j, double h) {
        return static_cast<double>(j) * h;
    }
    double y(int i, double h) {
        return 1.0 - static_cast<double>(i) * h;
    }


    double f(double x, double y) {
        return 0.0-2.0 * (x * (x - 1.0) + y * (y - 1.0));
    }

	vector<vector<double>> fillingVectorK(vector<vector<double>>& matrix, int N) {
		vector<vector<double>> temp((N-1)*(N-1));
        int kvB = 0;
        int kvE = N - 1;
        for (int i = 0; i < (N - 1) * (N - 1); i++)
        {
            for (int j = 0; j < (N - 1) * (N - 1); j++) {
                if (i + 1 > kvE)
                {
                    kvB += N - 1;
                    kvE += N - 1;
                }
                if (j < kvE && j >= kvB) {
                    if (i == j) temp[i].push_back(4.0);
                    else {
                        if (abs(i - j) == 1) temp[i].push_back(-1.0);
                        else temp[i].push_back(0.0);
                    }
                }
                else
                {
                    if (abs(i - j) == N - 1) temp[i].push_back(-1.0);
                    else temp[i].push_back(0.0);
                }
                
                /*if (i == j) {
                    temp[i].push_back(4.0);
                }
                else if ((abs(j - i) == 1 or abs(j - i) == N-1)) temp[i].push_back(-1.0);
                else
                {
                    temp[i].push_back(0);
                }*/
            }
        }
        /*for (int k = N-2; k < (N - 1) * (N - 1); k += (N - 1))
        {
            if (k + (N - 1) >= (N - 1) * (N - 1)) break;
            temp[k][k + 1] = 0.0;
            temp[k + 1][k] = 0.0;
        }*/
        return temp;
	}
    
    vector<double> fillVectorX(vector<double>& vec, int N) {
        vector<vector<double>> temp;
        
        
        for (int i = 0; i < N+1; i++)
        {
            vector<double> temp_2;
            for (int j = 0; j < N+1; j++)
            {
                temp_2.push_back(f(x(i, h), y(j, h)));
            }
            temp.push_back(temp_2);
        }
        

        vector<double> temp_3;
        for (int i = 1; i < N; i++)
        {
            for (int j = 1; j < N; j++)
            {
                temp_3.push_back((h * h / 12.0) * (6.0 * temp[i][j] + (temp[i][j - 1] + temp[i][j + 1] + temp[i + 1][j] + temp[i + 1][j - 1] + temp[i - 1][j] + temp[i - 1][j + 1])));
            }
        }
        return temp_3;
    }
}