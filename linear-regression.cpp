#include <stdio.h>
#include <iostream>
#include <cassert>
#include "linear-regression.hpp"
#include "utils.hpp"
#include "vector-ops.hpp"

// constructors and destructor
LinearRegression::LinearRegression() {}
LinearRegression::LinearRegression(std::string name) : name(name) {}
LinearRegression::~LinearRegression() {}

// public methods
LinearRegression &LinearRegression::fit(std::vector<std::vector<double>> X, std::vector<std::vector<double>> y, double lr, int epochs)
{
    /* all input arrays must be two dimensional */
    // validate input vectors
    assert((X.size() == y.size()) && (X.size() > 0) && (X[0].size() > 0) && (y[0].size() > 0));

    // store input dimensions
    int m = X.size(), n = X[0].size();

    // initialize zero weights - n x 1
    w = std::vector<std::vector<double>>(n, std::vector<double>(1, 0));

    // initialize bias to zero
    b = 0;

    // initialize used variables
    std::vector<std::vector<double>> y_hat, error, dw;
    double db;

    for (int epoch = 0; epoch < epochs; epoch++)
    {
        y_hat = matmul(X, w) + b;

        // std::cout << "y_hat\n";
        // print2d(y_hat);

        error = y_hat - y;

        // std::cout << "error\n";
        // print2d(error);

        dw = matmul(transpose(X), error) / m;

        // std::cout << "dw\n";
        // print2d(dw);

        db = sum(error) / m;
        w = w - lr * dw;

        // std::cout << "new w\n";
        // print2d(w);

        b = b - lr * db;
    }

    return *this;
}