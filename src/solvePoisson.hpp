/*----------------------------------------------------------------------------*\
  ME5311, Spring 2020
  Term Project
  Julian M. Toumey
--------------------------------------------------------------------------------
File
    solvePoisson.hpp

Description
    Solves the pressure Poisson equation for a given input field (typically 
    divergence of velocity). 
\*----------------------------------------------------------------------------*/

#ifndef SOLVEPOISSON_H
#define SOLVEPOISSON_H

#include "array2d.hpp"
#include "transforms.hpp"

Array2d solvePoisson(Array2d& inputArray, double dx, double dy)
{
    int M = inputArray.rowSize();
    int N = inputArray.colSize();

    Array2d D(M, N, 0.0);
 
    // Modify Transforms to return an Array2d instance
    Transforms::dct2d(inputArray, D);

    Array2d p(M, N, 0.0);

    for (int ii = 0; ii < M; ++ii)
    {
        for (int jj = 0; jj < N; ++jj)
        {
            p(ii, jj) = D(ii, jj)/((2*std::cos(M_PI*ii/M)-2)/(dx*dx) + (2*std::cos(M_PI*jj/N)-2)/(dy*dy));
        }
    }

    p(0, 0) = 0.0;
    Array2d pressure(M, N, 0.0);
    Transforms::idct2d(p, pressure);

    return (pressure);
}

#endif // SOLVEPOISSON_H