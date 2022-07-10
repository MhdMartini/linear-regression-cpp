
#include <iostream>
#include "utils.hpp"

std::vector<std::vector<double>> getUniformData(int m, int n)
{
    /*get m x n features of uniform distribution*/
    std::vector<std::vector<double>> X(m, std::vector<double>(n));
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
            X[i][j] = ((double)rand() / (RAND_MAX));
    }
    return X;
}
