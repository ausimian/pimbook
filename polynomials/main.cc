#include <stdio.h>
#include "polynomial.h"

int main(int, char**)
{
    polynomial f1({ 1.0, -3.0, 7.0 });

    points pts;
    for (auto x : {-2.0, 4.0, 10.0})
        pts.push_back({ x, f1(x) });

    polynomial f2(pts);

    for (auto c : f2.coefficients())
        printf("%f\n", c);
}