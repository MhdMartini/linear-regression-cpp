// #ifndef UTILS_H
// #define UTILS_H
#include <vector>
#include <iostream>

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

template <typename T>
void print1d(std::vector<T> v, std::string end = "\n")
{
    std::cout << "[";
    for (int i = 0; i < v.size(); i++)
    {
        std::cout << v[i];
        if (i == v.size() - 1)
            std::cout << "]" + end;
        else
            std::cout << ", ";
    }
}

template <typename T>
void print2d(std::vector<std::vector<T>> v)
{
    std::cout << "[";
    for (int i = 0; i < v.size(); i++)
    {
        print1d(v[i], "");
        if (!(i == v.size() - 1))
            std::cout << std::endl;
    }
    std::cout << "]\n";
}
// #endif // UTILS_H