#pragma once

#include <vector>

/*--
    Coefficients are expressed in order of increasing powers of x.
--*/
using coeffs = std::vector<double>;
using points = std::vector<std::pair<double,double>>;

class polynomial
{
public:
    polynomial() 
    {}

    explicit polynomial(coeffs cfs)
        : mCoeffs(std::move(cfs))
    {
    }

    explicit polynomial(std::initializer_list<double> cfs)
        : mCoeffs(cfs)
    {
    }

    explicit polynomial(const points& pts)
    {
        polynomial sum;
        for (size_t i = 0; i < pts.size(); ++i)
        {
            auto [xi, yi] = pts[i];

            polynomial term { 1.0 };
            for (size_t j = 0; j < pts.size(); ++j)
            {
                if (j == i)
                    continue;

                auto xj = pts[j].first;
                term *= polynomial { -xj / (xi - xj) , 1.0 / (xi - xj) };
            }

            sum += term * polynomial { yi };
        }

        mCoeffs = std::move(sum.mCoeffs);
    }

    double operator()(double x) const
    {
        double result = 0;

        double powx = 1.0;
        for (auto cf : mCoeffs)
        {
            result += cf * powx;
            powx *= x;
        }

        return result;
    };

    polynomial operator+(const polynomial& rhs) const
    {
        return polynomial(zip(mCoeffs, rhs.mCoeffs, [](double c1, double c2) { return c1 + c2;}));
    } 

    polynomial operator-(const polynomial& rhs) const
    {
        return polynomial(zip(mCoeffs, rhs.mCoeffs, [](double c1, double c2) { return c1 - c2;}));
    }

    polynomial operator*(const polynomial& rhs) const
    {
        return polynomial(multiply(mCoeffs, rhs.mCoeffs));
    }

    polynomial& operator+=(const polynomial& rhs)
    {
        mCoeffs = zip(mCoeffs, rhs.mCoeffs, [](double c1, double c2) { return c1 + c2;});
        return *this;
    }

    polynomial& operator*=(const polynomial& rhs)
    {
        mCoeffs = multiply(mCoeffs, rhs.mCoeffs);
        return *this;
    }

    int degree() const { return degree(mCoeffs); }

    const coeffs& coefficients() const { return mCoeffs; }

  private:
    template<typename F> static coeffs zip(const coeffs &lhs, const coeffs &rhs, F&& f)
    {
        coeffs cfs(std::max(lhs.size(), rhs.size()));

        auto itRes = cfs.begin();
        auto itLhs = lhs.begin();
        auto itRhs = rhs.begin();
        while (itLhs != lhs.end() || itRhs != rhs.end())
        {
            if (itLhs != lhs.end())
                *itRes = *itLhs++;
            if (itRhs != rhs.end())
                *itRes = f(*itRes, *itRhs++);
            ++itRes; 
        }

        return cfs;
    }

    static coeffs multiply(const coeffs &lhs, const coeffs &rhs)
    {
        coeffs res(degree(lhs) + degree(rhs) + 1);

        for (size_t i = 0; i < lhs.size(); ++i)
            for (size_t j = 0; j < rhs.size(); ++j)
                res[i + j] += lhs[i] * rhs[j];

        return res;
    }

    static int degree(const coeffs& cfs)
    {
        return static_cast<int>(cfs.size()) - 1;
    }

    coeffs mCoeffs;
};