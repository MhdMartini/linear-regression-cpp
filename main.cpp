#include <string>
#include "utils.hpp"
#include "vector-ops.hpp"
#include "linear-regression.hpp"

// #include "utils.hpp"

int main()
{

    // create features
    int m = 1000, n = 3;
    std::vector<std::vector<double>> X = getUniformData(m, n);

    // create random weights and bias
    std::vector<std::vector<double>> w = getUniformData(n, 1);
    double b = getUniformData(1, 1)[0][0];

    // print weights and bias
    std::cout << "weights:\n";
    print2d(w);
    std::cout << "bias:\n"
              << b << std::endl;

    // calculate y
    std::vector<std::vector<double>> y = matmul(X, w) + b;

    // run linear regression
    LinearRegression lr;
    lr.fit(X, y, 0.01, 9000);

    std::cout << "\nleared weights:\n";
    print2d(lr.w);

    std::cout << "leared bias:\n";
    std::cout << b << std::endl;

    return 0;
}