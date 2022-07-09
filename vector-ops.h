// #ifndef VECTOROPS_H
// #define VECTOROPS_H

#include <vector>
#include <stdexcept>
#include <iostream>

/* scaler and vector */
// scalar + 1d vector
template <typename L, typename R>
std::vector<double> operator+(L lhs, const std::vector<R> &rhs)
{
    std::vector<double> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = rhs[i] + lhs;
    return out;
}

// 1d vector + scalar
template <typename L, typename R>
std::vector<double> operator+(const std::vector<L> &lhs, R rhs) { return rhs + lhs; }

// 1d vector += scalar
template <typename L, typename R>
std::vector<double> operator+=(std::vector<L> &lhs, R rhs)
{
    for (int i = 0; i < lhs.size(); i++)
        lhs[i] += rhs;
    return lhs;
}

// scalar - 1d vector
template <typename L, typename R>
std::vector<double> operator-(L lhs, const std::vector<R> &rhs)
{
    std::vector<double> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = lhs - rhs[i];
    return out;
}

// 1d vector - scalar
template <typename L, typename R>
std::vector<double> operator-(const std::vector<L> &lhs, R rhs)
{
    std::vector<double> out(lhs.size());
    for (int i = 0; i < lhs.size(); i++)
        out[i] = lhs[i] - rhs;
    return out;
}

// 1d vector -= scalar
template <typename L, typename R>
std::vector<double> operator-=(std::vector<L> &lhs, R rhs)
{
    for (int i = 0; i < lhs.size(); i++)
        lhs[i] -= rhs;
    return lhs;
}

// scalar * 1d vector
template <typename L, typename R>
std::vector<double> operator*(L lhs, const std::vector<R> &rhs)
{
    std::vector<double> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = rhs[i] * lhs;
    return out;
}

// 1d vector * scalar
template <typename L, typename R>
std::vector<double> operator*(const std::vector<L> &lhs, R rhs) { return rhs * lhs; }

// 1d vector *= scalar
template <typename L, typename R>
std::vector<double> operator*=(std::vector<L> &lhs, R rhs)
{
    for (int i = 0; i < lhs.size(); i++)
        lhs[i] *= rhs;
    return lhs;
}

// scalar / 1d vector
template <typename L, typename R>
std::vector<double> operator/(L lhs, const std::vector<R> &rhs)
{
    std::vector<double> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = lhs / rhs[i];
    return out;
}

// 1d vector / scalar
template <typename L, typename R>
std::vector<double> operator/(const std::vector<L> &lhs, R rhs)
{
    std::vector<double> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = rhs[i] / lhs;
    return out;
}

// 1d vector /= scalar
template <typename L, typename R>
std::vector<double> operator/=(std::vector<L> &lhs, R rhs)
{
    for (int i = 0; i < lhs.size(); i++)
        lhs[i] /= rhs;
    return lhs;
}

/* 2d vector and scalar */
// 2d vector + scalar
template <typename L, typename R>
std::vector<std::vector<double>> operator+(const std::vector<std::vector<L>> &lhs, R rhs)
{
    std::vector<std::vector<double>> out(lhs.size());
    for (int i = 0; i < lhs.size(); i++)
        out[i] = lhs[i] + rhs;
    return out;
}

// scalar + 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator+(L lhs, const std::vector<std::vector<R>> &rhs)
{
    std::vector<std::vector<double>> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = lhs + rhs[i];
    return out;
}

// 2d vector += scalar
template <typename L, typename R>
std::vector<std::vector<double>> operator+=(std::vector<std::vector<L>> &lhs, R rhs)
{
    for (int i = 0; i < lhs.size(); i++)
        lhs[i] += rhs;
    return lhs;
}

// 2d vector - scalar
template <typename L, typename R>
std::vector<std::vector<double>> operator-(const std::vector<std::vector<L>> &lhs, R rhs)
{
    std::vector<std::vector<double>> out(lhs.size());
    for (int i = 0; i < lhs.size(); i++)
        out[i] = lhs[i] - rhs;
    return out;
}

// scalar - 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator-(L lhs, const std::vector<std::vector<R>> &rhs)
{
    std::vector<std::vector<double>> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = lhs - rhs[i];
    return out;
}

// 2d vector -= scalar
template <typename L, typename R>
std::vector<std::vector<double>> operator-=(std::vector<std::vector<L>> &lhs, R rhs)
{
    for (int i = 0; i < lhs.size(); i++)
        lhs[i] -= rhs;
    return lhs;
}

// 2d vector * scalar
template <typename L, typename R>
std::vector<std::vector<double>> operator*(const std::vector<std::vector<L>> &lhs, R rhs)
{
    std::vector<std::vector<double>> out(lhs.size());
    for (int i = 0; i < lhs.size(); i++)
        out[i] = lhs[i] * rhs;
    return out;
}

// scalar * 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator*(L lhs, const std::vector<std::vector<R>> &rhs)
{
    std::vector<std::vector<double>> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = rhs[i] * lhs;
    return out;
}

// 2d vector *= scalar
template <typename L, typename R>
std::vector<std::vector<double>> operator*=(std::vector<std::vector<L>> &lhs, R rhs)
{
    for (int i = 0; i < lhs.size(); i++)
        lhs[i] *= rhs;
    return lhs;
}

// 2d vector / scalar
template <typename L, typename R>
std::vector<std::vector<double>> operator/(const std::vector<std::vector<L>> &lhs, R rhs)
{
    std::vector<std::vector<double>> out(lhs.size());
    for (int i = 0; i < lhs.size(); i++)
        out[i] = lhs[i] / rhs;
    return out;
}

// scalar / 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator/(L lhs, const std::vector<std::vector<R>> &rhs)
{
    std::vector<std::vector<double>> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = lhs / rhs[i];
    return out;
}

// 2d vector /= scalar
template <typename L, typename R>
std::vector<std::vector<double>> operator/=(std::vector<std::vector<L>> &lhs, R rhs)
{
    for (int i = 0; i < lhs.size(); i++)
        lhs[i] /= rhs;
    return lhs;
}

/* 1d vector and 1d vector */
// 1d vector + 1d vector
template <typename L, typename R>
std::vector<double> operator+(const std::vector<L> &lhs, const std::vector<R> &rhs)
{
    // make sure both vectors have the same size
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("vector sizes do not match");

    std::vector<double> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = rhs[i] + lhs[i];
    return out;
}

// 1d vector += 1d vector
template <typename L, typename R>
std::vector<double> operator+=(std::vector<L> &lhs, const std::vector<R> &rhs)
{
    // make sure both vectors have the same size
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("vector sizes do not match");

    for (int i = 0; i < rhs.size(); i++)
        lhs[i] += rhs[i];
    return lhs;
}

// 1d vector - 1d vector
template <typename L, typename R>
std::vector<double> operator-(const std::vector<L> &lhs, const std::vector<R> &rhs)
{
    // make sure both vectors have the same size
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("vector sizes do not match");

    std::vector<double> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = rhs[i] - lhs[i];
    return out;
}

// 1d vector -= 1d vector
template <typename L, typename R>
std::vector<double> operator-=(std::vector<L> &lhs, const std::vector<R> &rhs)
{
    // make sure both vectors have the same size
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("vector sizes do not match");

    for (int i = 0; i < rhs.size(); i++)
        lhs[i] -= rhs[i];
    return lhs;
}

// 1d vector * 1d vector
template <typename L, typename R>
std::vector<double> operator*(const std::vector<L> &lhs, const std::vector<R> &rhs)
{
    // make sure both vectors have the same size
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("vector sizes do not match");

    std::vector<double> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = rhs[i] * lhs[i];
    return out;
}

// 1d vector *= 1d vector
template <typename L, typename R>
std::vector<double> operator*=(std::vector<L> &lhs, const std::vector<R> &rhs)
{
    // make sure both vectors have the same size
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("vector sizes do not match");

    for (int i = 0; i < rhs.size(); i++)
        lhs[i] *= rhs[i];
    return lhs;
}

// 1d vector / 1d vector
template <typename L, typename R>
std::vector<double> operator/(const std::vector<L> &lhs, const std::vector<R> &rhs)
{
    // make sure both vectors have the same size
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("vector sizes do not match");

    std::vector<double> out(rhs.size());
    for (int i = 0; i < rhs.size(); i++)
        out[i] = lhs[i] / rhs[i];
    return out;
}

// 1d vector /= 1d vector
template <typename L, typename R>
std::vector<double> operator/=(std::vector<L> &lhs, const std::vector<R> &rhs)
{
    // make sure both vectors have the same size
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("vector sizes do not match");

    for (int i = 0; i < rhs.size(); i++)
        lhs[i] /= rhs[i];
    return lhs;
}

// 2d vector + 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator+(const std::vector<std::vector<L>> &lhs, const std::vector<std::vector<R>> &rhs)
{
    // make sure both vectors have the same size
    if ((lhs.size() != rhs.size()) | (lhs[0].size() != rhs[0].size()))
        throw std::invalid_argument("vector sizes do not match");

    std::vector<std::vector<double>> out = lhs;
    for (int i = 0; i < lhs.size(); i++)
        out[i] += rhs[i];
    return out;
}

// 2d vector += 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator+=(std::vector<std::vector<L>> &lhs, const std::vector<std::vector<R>> &rhs)
{
    // make sure both vectors have the same size
    if ((lhs.size() != rhs.size()) | (lhs[0].size() != rhs[0].size()))
        throw std::invalid_argument("vector sizes do not match");

    for (int i = 0; i < lhs.size(); i++)
        lhs[i] += rhs[i];
    return lhs;
}

// 2d vector - 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator-(const std::vector<std::vector<L>> &lhs, const std::vector<std::vector<R>> &rhs)
{
    // make sure both vectors have the same size
    if ((lhs.size() != rhs.size()) | (lhs[0].size() != rhs[0].size()))
        throw std::invalid_argument("vector sizes do not match");

    std::vector<std::vector<double>> out = lhs;
    for (int i = 0; i < lhs.size(); i++)
        out[i] -= rhs[i];
    return out;
}

// 2d vector -= 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator-=(std::vector<std::vector<L>> &lhs, const std::vector<std::vector<R>> &rhs)
{
    // make sure both vectors have the same size
    if ((lhs.size() != rhs.size()) | (lhs[0].size() != rhs[0].size()))
        throw std::invalid_argument("vector sizes do not match");

    for (int i = 0; i < lhs.size(); i++)
        lhs[i] -= rhs[i];
    return lhs;
}

// 2d vector * 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator*(const std::vector<std::vector<L>> &lhs, const std::vector<std::vector<R>> &rhs)
{
    // make sure both vectors have the same size
    if ((lhs.size() != rhs.size()) | (lhs[0].size() != rhs[0].size()))
        throw std::invalid_argument("vector sizes do not match");

    std::vector<std::vector<double>> out(lhs.size(), std::vector<double>(lhs[0].size()));
    for (int i = 0; i < lhs.size(); i++)
        for (int j = 0; j < lhs[0].size(); j++)
            out[i][j] = lhs[i][j] * rhs[i][j];
    return out;
}

// 2d vector *= 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator*=(std::vector<std::vector<L>> &lhs, const std::vector<std::vector<R>> &rhs)
{
    // make sure both vectors have the same size
    if ((lhs.size() != rhs.size()) | (lhs[0].size() != rhs[0].size()))
        throw std::invalid_argument("vector sizes do not match");

    for (int i = 0; i < lhs.size(); i++)
        for (int j = 0; j < lhs[0].size(); j++)
            lhs[i][j] *= rhs[i][j];
    return lhs;
}

// 2d vector / 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator/(const std::vector<std::vector<L>> &lhs, const std::vector<std::vector<R>> &rhs)
{
    // make sure both vectors have the same size
    if ((lhs.size() != rhs.size()) | (lhs[0].size() != rhs[0].size()))
        throw std::invalid_argument("vector sizes do not match");

    std::vector<std::vector<double>> out(lhs.size(), std::vector<double>(lhs[0].size()));
    for (int i = 0; i < lhs.size(); i++)
        for (int j = 0; j < lhs[0].size(); j++)
            out[i][j] = lhs[i][j] / rhs[i][j];
    return out;
}

// 2d vector /= 2d vector
template <typename L, typename R>
std::vector<std::vector<double>> operator/=(std::vector<std::vector<L>> &lhs, const std::vector<std::vector<R>> &rhs)
{
    // make sure both vectors have the same size
    if ((lhs.size() != rhs.size()) | (lhs[0].size() != rhs[0].size()))
        throw std::invalid_argument("vector sizes do not match");

    for (int i = 0; i < lhs.size(); i++)
        for (int j = 0; j < lhs[0].size(); j++)
            lhs[i][j] /= rhs[i][j];
    return lhs;
}

// #endif // VECTOROPS_H