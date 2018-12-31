#pragma once

#include <array>

template <typename T, size_t r, size_t c = r>
using matrix = std::array<std::array<T, c>, r>;

template <typename T, size_t rows, size_t cols, typename = std::enable_if_t<cols == rows + 1>>
bool to_upper_tri(matrix<T, rows, cols>& m)
{
    const auto swap_row = [&](size_t i, size_t j) 
    { 
        for (size_t k = 0; k < cols; ++k)
            std::swap(m[i][k], m[j][k]);
    };

    for (size_t n = 0; n < rows; ++n)
    {
        size_t i = n;
        while (i < rows && !m[i][n])
            ++i;
        if (i == rows)
            return false;
        if (i > n)
            swap_row(i, n);

        const auto &r1 = m[n];
        for (i = n + 1; i < rows; ++i)
        {
            auto& r2 = m[i];
            if (!r2[n])
                continue;

            const auto lcm = std::lcm(r1[n],r2[n]);
            const auto r1f = lcm / r1[n];
            const auto r2f = lcm / r2[n];

            for (size_t j = n; j < cols; ++j)
                r2[j] = r2[j] * r2f - r1[j] * r1f;
        }
    }

    return true;
}

template <typename T, size_t rows, size_t cols, typename = std::enable_if_t<cols == rows + 1>>
bool solve(const matrix<T, rows, cols>& input, std::array<T, rows>& result)
{
    auto utri = input;
    if (!to_upper_tri(utri))
        return false;

    result[rows - 1] = utri[rows - 1][rows] / utri[rows - 1][rows - 1];
    
    for (size_t j = rows - 2; j < rows; --j)
    {
        auto rhs = utri[j][rows];
        for (size_t k = rows - 1; k > j; --k)
            rhs -= result[k] * utri[j][k];
        result[j] = rhs / utri[j][j];
    }

    return true;
}
