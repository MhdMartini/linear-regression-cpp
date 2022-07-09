#include <stdio.h>
#include <iostream>
#include "linear-regression.h"

// constructors and destructor
LinearRegression::LinearRegression() {}
LinearRegression::LinearRegression(std::string name) : name(name) {}
LinearRegression::~LinearRegression() {}

// public methods
LinearRegression &LinearRegression::fit(std::vector<std::vector<double>> X, std::vector<double> y, int epochs)
{
    int m = X.size();
    int n = X[0].size();

    // initialize weights
    std::vector<double> w(n);

    return *this;
}