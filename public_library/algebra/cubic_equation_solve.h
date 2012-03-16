//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// ditto/public_library/algebra/cubic_equation_solve.h
// Copyright 2012, Chenfanfu Jiang
// cubic equation solver using Cardano's method
// A solver for ax^3+bx^2+c^x+d=0
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#ifndef DITTO_PUBLIC_LIBRARY_ALGEBRA_CUBIC_EQUATION_SOLVE_H
#define DITTO_PUBLIC_LIBRARY_ALGEBRA_CUBIC_EQUATION_SOLVE_H

#include <cmath>
#include <iostream>

namespace ditto { namespace algebra {

#define CUBIC_SOLVER_THIRD 0.333333333333333
#define CUBIC_SOLVER_ROOTTHREE 1.73205080756888        

// this function returns the cube root if x were a negative number aswell
template<class T>
T cube_root(T x)
{
    if (x < 0)
        return -std::pow(-x, CUBIC_SOLVER_THIRD);
    else
        return std::pow(x, CUBIC_SOLVER_THIRD);
}

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
// Function cubic_equation_solve
// Return value: 1: three equal real roots
//               2: three (distinct?) real roots
//               3: 1 real root and 2 complex roots
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
template<class T, class ComplexNumberType>
int cubic_equation_solve(T a, T b, T c, T d, ComplexNumberType &first_root, ComplexNumberType &second_root, ComplexNumberType &third_root)
{
    int solution_nature;

    // find the discriminant
    T f, g, h;
    f = (3 * c / a - std::pow(b, 2) / std::pow(a, 2)) / 3;
    g = (2 * std::pow(b, 3) / std::pow(a, 3) - 9 * b * c / std::pow(a, 2) + 27 * d / a) / 27;
    h = std::pow(g, 2) / 4 + std::pow(f, 3) / 27;
    // evaluate discriminant
    if (f == 0 && g == 0 && h == 0)
    {
        // 3 equal roots
        T x;
        // when f, g, and h all equal 0 the roots can be found by the following line
        x = -cube_root(d / a);

        solution_nature = 1;
        first_root(0) = x;
        second_root(0) = x;
        third_root(0) = x;
    }
    else if (h <= 0)
    {
        // 3 real roots
        T q, i, j, k, l, m, n, p;
        // complicated maths making use of the method
        i = std::pow(std::pow(g, 2) / 4 - h, 0.5);
        j = cube_root(i);
        k = acos(-(g / (2 * i)));
        m = cos(k / 3);
        n = CUBIC_SOLVER_ROOTTHREE * sin(k / 3);
        p = -(b / (3 * a));

        solution_nature = 2;
        first_root(0) = 2 * j * m + p;
        second_root(0) = -j * (m + n) + p;
        third_root(0) = -j * (m - n) + p;
    }
    else if (h > 0)
    {
        // 1 real root and 2 complex roots
        T r, s, t, u, p;
        // complicated maths making use of the method
        r = -(g / 2) + std::pow(h, 0.5);
        s = cube_root(r);
        t = -(g / 2) - std::pow(h, 0.5);
        u = cube_root(t);
        p = -(b / (3 * a));

        solution_nature = 3;
        first_root(0) = (s + u) + p;
        second_root(0) =  -(s + u) / 2 + p; second_root(1) =  (s - u) * CUBIC_SOLVER_ROOTTHREE / 2;
        third_root(0) =  -(s + u) / 2 + p; third_root(1) =  -(s - u) * CUBIC_SOLVER_ROOTTHREE / 2;
    }

    return solution_nature;
}


}} // end of namespaces
#endif 
