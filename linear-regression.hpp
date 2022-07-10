// #ifndef LINEARREGRESSION_H
// #define LINEARREGRESSION_H

#include <vector>
#include <string>

class LinearRegression
{
public:
    std::string name = "LinearRegression";
    std::vector<std::vector<double>> w;
    double b;

public:
    LinearRegression();
    LinearRegression(std::string name);
    ~LinearRegression();

public:
    LinearRegression &fit(std::vector<std::vector<double>> X, std::vector<std::vector<double>> y, double lr = 0.01, int epochs = 1000);
};

// #endif // LINEARREGRESSION_H