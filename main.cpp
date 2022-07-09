#include <string>
#include "linear-regression.h"
#include "utils.h"
#include "vector-ops.h"

int main()
{
    LinearRegression lr;
    int m = 2, n = 4;

    std::vector<std::vector<double>> X = getUniformData(m, n);
    print2d(X);
    print2d(2 / X);

    return 0;
}