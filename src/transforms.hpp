/*----------------------------------------------------------------------------*\
  ME5311, Spring 2020
  Term Project
  Julian M. Toumey
--------------------------------------------------------------------------------
File
    transforms.hpp

Description
    Collection of functions for computing 2D forward and reverse discrete cosine 
    transforms. 

    Note: The COEFFS macro in the iDCT2D was borrowed directly from 
    E. Mikulic's article (https://unix4lyfe.org/dct/). Large parts of the iDCT2D 
    in this file are taken from the listings on the above page.
\*----------------------------------------------------------------------------*/

#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include "array2d.hpp"
#include <cmath>

namespace Transforms
{

double dct(double input, double length, int ii, int jj)
{
    double output = input*std::cos((M_PI/length)*(jj + 0.5)*ii);
    return output;
}

void dct2d(Array2d& inputSignal, Array2d& outputSignal)
{
    int numPointsk1 = inputSignal.rowSize();
    int numPointsk2 = inputSignal.colSize();

    Array2d tempDct(numPointsk1, numPointsk2, 0.0);

    const double eps = 1e-12;
    double alphak1;
    double alphak2;

    // Build look-up table
    double* cosResRow = new double[numPointsk2*numPointsk2];
    for (int jj = 0; jj < numPointsk2; ++jj)
    {
        for (int kk = 0; kk < numPointsk2; ++kk)
        {
            cosResRow[jj*numPointsk2 + kk] = std::cos((M_PI/numPointsk2)*(kk + 0.5)*jj);
        }
    }

    double* cosResCol = new double[numPointsk1*numPointsk1];
    for (int ii = 0; ii < numPointsk1; ++ii)
    {
        for (int kk = 0; kk < numPointsk1; ++kk)
        {
            cosResCol[ii*numPointsk1 + kk] = std::cos((M_PI/numPointsk1)*(kk + 0.5)*ii);
        }
    }


    // Peform a 1D DCT on every row
    for (int ii = 0; ii < numPointsk1; ++ii)
    {
        for (int jj = 0; jj < numPointsk2; ++jj)
        {
            for (int kk = 0; kk < numPointsk2; ++kk)
            {
                tempDct(ii, jj) += inputSignal(ii, kk)*cosResRow[jj*numPointsk2 + kk];
                // tempDct(ii, jj) += inputSignal(ii, kk)*std::cos((M_PI/numPointsk2)*(kk + 0.5)*jj);
            }
            if (jj == 0)
            {
                alphak1 = std::sqrt(1.0/(double) numPointsk2);
            }
            else 
            {
                alphak1 = std::sqrt(2.0/(double) numPointsk2);
            }
            tempDct(ii, jj) = tempDct(ii, jj)*alphak1;
        }
    }

    // Perform a 1D DCT on every column
    for (int jj = 0; jj < numPointsk2; ++jj)
    {
        for (int ii = 0; ii < numPointsk1; ++ii)
        {
            for (int kk = 0; kk < numPointsk1; ++kk)
            {
                outputSignal(ii, jj) += tempDct(kk, jj)*cosResCol[ii*numPointsk1 + kk];
                // outputSignal(ii, jj) += tempDct(kk, jj)*std::cos((M_PI/numPointsk1)*(kk + 0.5)*ii);
            }
            if (ii == 0)
            {
                alphak2 = std::sqrt(1.0/(double) numPointsk1);
            }
            else
            {
                alphak2 = std::sqrt(2.0/(double) numPointsk1);
            }
            outputSignal(ii, jj) = outputSignal(ii, jj)*alphak2;

            // Eliminate very small (essentially zero) values based on machine 
            // precision
            if (std::abs(outputSignal(ii, jj)) < eps)
            {
                outputSignal(ii, jj) = 0.0;
            }
        }
    }

    delete [] cosResRow;
    delete [] cosResCol;
}


void idct2d(Array2d& inputSignal, Array2d& outputSignal)
{
    int numPointsk1 = inputSignal.rowSize();
    int numPointsk2 = inputSignal.colSize();

    double* cosResRow = new double[numPointsk2*numPointsk2];
    for (int yy = 0; yy < numPointsk2; ++yy)
    {
        for (int vv = 0; vv < numPointsk2; ++vv)
        {
            cosResRow[yy*numPointsk2 + vv] = std::cos((M_PI*(2.0*(double)yy + 1)*
                    (double)vv)/(2.0*(double) numPointsk2));
            // std::cos((M_PI/numPointsk2)*(kk + 0.5)*jj);
        }
    }

    double* cosResCol = new double[numPointsk1*numPointsk1];
    for (int xx = 0; xx < numPointsk1; ++xx)
    {
        for (int uu = 0; uu < numPointsk1; ++uu)
        {
            cosResCol[xx*numPointsk1 + uu] = std::cos((M_PI*(2.0*(double)xx + 1)*(double)uu)/(2.0*(double) numPointsk1));
            // std::cos((M_PI/numPointsk1)*(kk + 0.5)*ii);
        }
    }

    
    #define COEFFS(Cu,Cv,uu,vv) { \
        if (uu == 0) Cu = 1.0 / std::sqrt((double) numPointsk1); else Cu = std::sqrt(2.0) / std::sqrt((double) numPointsk1); \
        if (vv == 0) Cv = 1.0 / std::sqrt((double) numPointsk2); else Cv = std::sqrt(2.0) / std::sqrt((double) numPointsk2); \
    }

    Array2d tempDct(numPointsk1, numPointsk2, 0.0);

    for (int xx = 0; xx < numPointsk1; ++xx)
    {
        for (int yy = 0; yy < numPointsk2; ++yy)
        {
            double zz = 0;

            for (int vv = 0; vv < numPointsk2; ++vv)
            {
                double S, q;
                double Cu, Cv;
                int uu = 1;
                COEFFS(Cu, Cv, uu, vv);

                S = inputSignal(xx, vv);
                // q = Cv*S*std::cos((M_PI*(2.0*(double)yy + 1)*
                //     (double)vv)/(2.0*(double) numPointsk2));
                q = Cv*S*cosResRow[yy*numPointsk2 + vv];

                zz = zz + q;

            }

            tempDct(xx, yy) = zz;
        }
    }

    for (int xx = 0; xx < numPointsk1; ++xx)
    {
        for (int yy = 0; yy < numPointsk2; ++yy)
        {
            double zz = 0;

            for (int uu = 0; uu < numPointsk1; ++uu)
            {
                double S, q;
                double Cu, Cv;
                int vv = 1;
                COEFFS(Cu, Cv, uu, vv);

                S = tempDct(uu, yy);
                q = Cu*S*cosResCol[xx*numPointsk1+uu];
                    // std::cos((M_PI*(2.0*(double)xx + 1)*(double)uu)/(2.0*(double) numPointsk1));

                zz = zz + q;

            }

            outputSignal(xx, yy) = zz;
        }
    }

    delete [] cosResRow;
    delete [] cosResCol;
}

} // End namespace Transforms

#endif // TRANSFORMS_H
