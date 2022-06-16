#pragma once
#include <vector>


namespace data {
	const int N = 10;
	const double l = 1;
	const double h = l / static_cast<double>(N);
    
    double x(int i, double h) {
        return static_cast<double>(i) * h;
    }
    double y(int j, double h) {
        return 1.0 - static_cast<double>(j) * h;
    }


    double f(double x, double y) {
        return -2.0 * (x * (x - 1) + y * (y - 1));
    }

	vector<vector<double>> fillingVectorK(int N) {
		vector<vector<double>> temp((N-1)*(N-1));
        for (int i = 0; i < (N - 1) * (N - 1); i++)
        {
            for (int j = 0; j < (N - 1) * (N - 1); j++) {
                if (i == j) {
                    temp[i].push_back(4.0);
                }
                else if ((abs(j - i) == 1 or abs(j - i) == 3)) temp[i].push_back(-1.0);
                else
                {
                    temp[i].push_back(0);
                }
            }
        }
        for (int i = 2; i < (N - 1) * (N - 1); i += (N - 1))
        {
            if (i + (N - 1) >= (N - 1) * (N - 1)) break;
            temp[i][i + 1] = 0.0;
            temp[i + 1][i] = 0.0;
        }
        return temp;
	}
    
    vector<double> fillVectorX(int N) {
        vector<vector<double>> temp((N + 1));
        vector<double> temp_2((N + 1), 0.0);
        temp[0] = temp_2;
        for (int i = 1; i < N; i++)
        {
            for (int j = 0; j < N+1; j++)
            {
                if (j == 0 or j == N) temp[i].push_back(0.0);
                else temp[i].push_back(f(x(i, h), y(j, h)));
            }
        }
        temp.back() = temp_2;

        vector<double> temp_3;
        for (int i = 1; i < N; i++)
        {
            for (int j = 1; j < N; j++)
            {
                temp_3.push_back((h * h / 12) * (6 * temp[i][j] + 2 * (temp[i][j - 1] + temp[i][j + 1] + temp[i + 1][j] + temp[i - 1][j])));
            }
        }
        return temp_3;
    }
}