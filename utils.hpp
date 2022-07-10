// #ifndef UTILS_H
// #define UTILS_H
#include <vector>
#include <string>
#include <iostream>

std::vector<std::vector<double>> getUniformData(int m, int n);

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

template <typename T>
std::vector<std::vector<T>> transpose(std::vector<std::vector<T>> v)
{
    int m = v.size(), n = v[0].size();
    std::vector<std::vector<T>> t(n, std::vector<T>(m));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
            t[i][j] = v[j][i];
    }
    return t;
}

template <typename T>
double sum(std::vector<T> &v)
// sum of all elements in 1d vector
{
    double s = 0;
    for (int i = 0; i < v.size(); i++)
        s += v[i];
    return s;
}

template <typename T>
double sum(std::vector<std::vector<T>> &v)
// sum of all elements in 2d vector
{
    double s = 0;
    for (int i = 0; i < v.size(); i++)
        s += sum(v[i]);
    return s;
}
// #endif // UTILS_H
