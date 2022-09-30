#pragma once
#include <vector>


namespace data {
	const int N = 5;
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
        int kvB = 0;
        int kvE = N - 1;
		vector<vector<double>> K(5);
        /*
        K[0] - last verh diag
        K[1] - first verh diag
        K[2] - central diag
        K[3] - first nizh diag
        K[5] - last nizh diag
        */
        K[3].push_back(0.0);
        K[4].insert(K[4].begin(), N - 1, 0.0);
        for (int i = 0; i < (N - 1) * (N - 1); i++)
        {
            for (int j = 0; j < (N - 1) * (N - 1); j++) {
                if (i + 1 > kvE)
                {
                    kvB += N - 1;
                    kvE += N - 1;
                }
                if (j < kvE && j >= kvB) {
                    if (i == j) K[2].push_back(4.0);
                    else {
                        if (j > i and abs(i - j) == 1) K[1].push_back(-1.0);
                        if (i > j and abs(i - j) == 1) K[3].push_back(-1.0);
                        //else temp[i].push_back(0.0);
                    }
                }
                else
                {
                    if (j > i and abs(i - j) == 1) K[1].push_back(0.0);
                    if (i > j and abs(i - j) == 1) K[3].push_back(0.0);
                    if (j > i and abs(i - j) == N - 1) K[0].push_back(-1.0);
                    if (i > j and abs(i - j) == N - 1) K[4].push_back(-1.0);
                    //else temp[i].push_back(0.0);
                }
            }
        }
        K[0].insert(K[0].end(), N - 1, 0.0);
        K[1].push_back(0.0);
        return K;
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