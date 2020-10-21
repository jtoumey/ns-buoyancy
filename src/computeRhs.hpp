/*----------------------------------------------------------------------------*\
  ME5311, Spring 2020
  Term Project
  Julian M. Toumey
--------------------------------------------------------------------------------
File
    computeRhs.hpp

Description
    Collection of functions for computing the discretized convection terms in 
    the Navier-Stokes equations. 
\*----------------------------------------------------------------------------*/

#ifndef COMPUTERHS_H
#define COMPUTERHS_H

#include "field.hpp"
#include "fieldCalculus.hpp"

Field computeRhsU(Field& uVel, Field& vVel, double dx, double dy)
{
    int M = vVel.rowSize();
    int N = uVel.colSize();
    
    int numGhostCells = uVel.ghostCellsPerBnd();    
    // Compute the uu flux
    Array2d fluxFu(M, N, 0.0, numGhostCells);
    for (int ii = uVel.startIdxX(); ii <= uVel.endIdxX() - 1; ++ii)
    {
        for (int jj = uVel.startIdxY(); jj <= uVel.endIdxY(); ++jj)
        {
            const double uCellCenter = (uVel(ii + 1, jj) + uVel(ii, jj))/2.0;
            fluxFu(ii, jj) = uCellCenter*uCellCenter;
        }
    }

    // Approximate the derivative of flux uu
    Array2d duuDx(M - 1, N, 0.0, numGhostCells);
    for (int ii = fluxFu.startIdxX(); ii <= fluxFu.endIdxX() - 1; ++ii)
    {
        for (int jj = fluxFu.startIdxY(); jj <= fluxFu.endIdxY(); ++jj)
        {
            duuDx(ii, jj) = (fluxFu(ii + 1, jj) - fluxFu(ii, jj))/dx;
        }
    }

    // Work arrays for cross flux calculation
    Array2d uCellCorner(M - 1, N, 0.0, numGhostCells);

    for (int ii = uVel.startIdxX() + 1; ii <= uVel.endIdxX() - 1; ++ii)
    {
        for (int jj = uVel.startIdxY(); jj <= uVel.endIdxY(); ++jj)
        {
            uCellCorner(ii - 1, jj) = (uVel(ii, jj) + uVel(ii, jj - 1))/2.0;
        }
    }

    Array2d vCellCorner(M - 1, N, 0.0, numGhostCells);
    // Again, we start from the first v velocity idx + 1 by convention for 
    // our linear interpolation
    for (int ii = vVel.startIdxX() + 1; ii <= vVel.endIdxX(); ++ii)
    {
        for (int jj = vVel.startIdxY(); jj <= vVel.endIdxY() - 1; ++jj)
        {
            vCellCorner(ii - 1, jj) = (vVel(ii, jj) + vVel(ii - 1, jj))/2.0;
        }
    }

    Array2d fluxGu(M - 1, N + 1, 0.0, numGhostCells);
    // Now, after playing that tedious game, we can compute the cross flux 
    // itself
    for (int ii = fluxGu.startIdxX(); ii <= fluxGu.endIdxX(); ++ii)
    {
        for (int jj = fluxGu.startIdxY(); jj <= fluxGu.endIdxY(); ++jj)
        {
            fluxGu(ii, jj) = uCellCorner(ii, jj)*vCellCorner(ii, jj);
        }
    }

    // Approximate the derivative of flux uv
    Array2d duvDy(M - 1, N, 0.0, numGhostCells);

    for (int ii = fluxGu.startIdxX(); ii <= fluxGu.endIdxX(); ++ii)
    {
        for (int jj = fluxGu.startIdxY(); jj <= fluxGu.endIdxY() - 1; ++jj)
        {
            duvDy(ii, jj) = (fluxGu(ii, jj + 1) - fluxGu(ii, jj))/dy;
        }
    }

    // Form the RHS
    // Note: Here, we allocate an rhs array that has M+1 x cells. Then, in the 
    // x-loop, the ii is modified by 1. In this way, we can pad out the rhsU 
    // array such that it matches the size of the uVel array, and after 
    // integration, we can add the new velocity.
    Field rhsU(M + 1, N, 0.0, numGhostCells);
    for (int ii = rhsU.startIdxX(); ii <= rhsU.endIdxX() - 1; ++ii)
    {
        for (int jj = rhsU.startIdxY(); jj <= rhsU.endIdxY(); ++jj)
        {
            rhsU(ii + 1, jj) = -duuDx(ii, jj) - duvDy(ii, jj);
        }
    }

    return (rhsU);
}

Field computeRhsV(Field& uVel, Field& vVel, double dx, double dy)
{
    int M = vVel.rowSize();
    int N = uVel.colSize();
    
    int numGhostCells = uVel.ghostCellsPerBnd();
    /// ** ---- **
    // Compute the vv flux
    Array2d fluxGv(M, N, 0.0, numGhostCells);
    for (int ii = vVel.startIdxX(); ii <= vVel.endIdxX(); ++ii)
    {
        for (int jj = vVel.startIdxY(); jj <= vVel.endIdxY() - 1; ++jj)
        {
            const double vCellCenter = (vVel(ii, jj + 1) + vVel(ii, jj))/2.0;
            fluxGv(ii, jj) = vCellCenter*vCellCenter;
        }
    }

    //Approximate the derivative of flux vv
    Array2d dvvDy(M, N - 1, 0.0, numGhostCells);
    for (int ii = fluxGv.startIdxX(); ii <= fluxGv.endIdxX(); ++ii)
    {
        for (int jj = fluxGv.startIdxY(); jj <= fluxGv.endIdxY() - 1; ++jj)
        {
            dvvDy(ii, jj) = (fluxGv(ii, jj + 1) - fluxGv(ii, jj))/dx;
        }
    }

    // Cross flux calculation
    Array2d uNwCorner(M, N - 1, 0.0, numGhostCells);
    for (int ii = uVel.startIdxX(); ii <= uVel.endIdxX(); ++ii)
    {
        for (int jj = uVel.startIdxY() + 1; jj <= uVel.endIdxY(); ++jj)
        {
            uNwCorner(ii, jj - 1) = (uVel(ii, jj) + uVel(ii, jj - 1))/2.0;
        }
    }

    Array2d vNwCorner(M, N - 1, 0.0, numGhostCells);
    for (int ii = vVel.startIdxX(); ii <= vVel.endIdxX(); ++ii)
    {
        for (int jj = vVel.startIdxY() + 1; jj <= vVel.endIdxY(); ++jj)
        {
            vNwCorner(ii, jj - 1) = (vVel(ii, jj) + vVel(ii - 1, jj))/2.0;
        }
    }

    Array2d fluxFv(M + 1, N - 1, 0.0, numGhostCells);
    for (int ii = fluxFv.startIdxX(); ii <= fluxFv.endIdxX(); ++ii)
    {
        for (int jj = fluxFv.startIdxY(); jj <= fluxFv.endIdxY(); ++jj)
        {
            fluxFv(ii, jj) = uNwCorner(ii, jj)*vNwCorner(ii, jj);
        }
    }

    // Approximate the derivative of flux vu 
    Array2d duvDx(M, N - 1, 0.0, numGhostCells);
    for (int ii = fluxFv.startIdxX(); ii <= fluxFv.endIdxX() - 1; ++ii)
    {
        for (int jj = fluxFv.startIdxY(); jj <= fluxFv.endIdxY(); ++jj)
        {
            duvDx(ii, jj) = (fluxFv(ii + 1, jj) - fluxFv(ii, jj))/dx;
        }
    }

    // Form the RHS
    // As with rhsU, pad out the array so it matches the size of uVel
    Field rhsV(M, N + 1, 0.0, numGhostCells);
    for (int ii = rhsV.startIdxX(); ii <= rhsV.endIdxX(); ++ii)
    {
        for (int jj = rhsV.startIdxY(); jj <= rhsV.endIdxY() - 1; ++jj)
        {
            rhsV(ii, jj + 1) = -duvDx(ii, jj) - dvvDy(ii, jj);
        }
    }
    return (rhsV);
}

///
Field computeRhsTheta(Field& theta, Field& uVel, Field& vVel, double dx, double dy)
{
    int M = theta.rowSize();
    int N = theta.colSize();
    
    int numGhostCells = (theta.totalCellsX() - M)/2;

    Array2d fluxFth(M, N, 0.0, numGhostCells);
    for (int ii = theta.startIdxX(); ii <= theta.endIdxX(); ++ii)
    {
        for (int jj = theta.startIdxY(); jj <= theta.endIdxY(); ++jj)
        {
            const double thetaWestface = (theta(ii - 1, jj) + theta(ii, jj))/2.0;

            fluxFth(ii, jj) = thetaWestface*uVel(ii, jj);
        }
    }

    Array2d fluxGth(M, N, 0.0, numGhostCells);
    for (int ii = theta.startIdxX(); ii <= theta.endIdxX(); ++ii)
    {
        for (int jj = theta.startIdxY(); jj <= theta.endIdxY(); ++jj)
        {
            const double thetaSouthface = (theta(ii, jj - 1) + 
                                           theta(ii, jj))/2.0;

            // Small, unfortunate, obfuscation: relying on ii, jj to pull 
            // correct velocity here. 
            fluxGth(ii, jj) = thetaSouthface*vVel(ii, jj);
        }
    }


    Array2d duThetaDx(M, N, 0.0, numGhostCells);
    Array2d dvThetaDy(M, N, 0.0, numGhostCells);

    for (int ii = duThetaDx.startIdxX(); ii <= duThetaDx.endIdxX(); ++ii)
    {
        for (int jj = duThetaDx.startIdxY(); jj <= duThetaDx.endIdxY(); ++jj)
        {
            duThetaDx(ii, jj) = (fluxFth(ii + 1, jj) - fluxFth(ii, jj))/dx;
        }
    }

    for (int ii = dvThetaDy.startIdxX(); ii <= dvThetaDy.endIdxX(); ++ii)
    {
        for (int jj = dvThetaDy.startIdxY(); jj <= dvThetaDy.endIdxY(); ++jj)
        {
            dvThetaDy(ii, jj) = (fluxGth(ii, jj + 1) - fluxGth(ii, jj))/dy;
        }
    }
    
    Field rhsTheta = duThetaDx*(-1) - dvThetaDy;

    return (rhsTheta);
}

Field computeRhsViscU(Field& uVel, double dx, double dy, double prandtlNumber)
{
    int M = uVel.rowSize();
    int N = uVel.colSize();
    int numGhostCells = uVel.ghostCellsPerBnd();
    Field rhsViscU(M, N, 0.0, numGhostCells);

    rhsViscU = (laplacianx(uVel, dx, dy))*prandtlNumber;

    return (rhsViscU);
}

Field computeRhsViscV(Field& vVel, double dx, double dy, double prandtlNumber)
{
    int M = vVel.rowSize();
    int N = vVel.colSize();
    int numGhostCells = vVel.ghostCellsPerBnd();
    Field rhsViscV(M, N, 0.0, numGhostCells);

    rhsViscV = (laplaciany(vVel, dx, dy))*prandtlNumber;

    return (rhsViscV);
}

Field computeRhsViscTheta(Field& theta, double dx, double dy, double prandtlNumber)
{
    int M = theta.rowSize();
    int N = theta.colSize();
    int numGhostCells = theta.ghostCellsPerBnd();
    Field rhsViscTheta(M, N, 0.0, numGhostCells);

    rhsViscTheta = (laplacianTheta(theta, dx, dy));

    return (rhsViscTheta);
}

Field computeRhsBuoyancyV(Field& theta, double rayleighNum, double prandtlNum)
{
    const double raTimesPr = rayleighNum*prandtlNum;
    int M = theta.rowSize();
    int N = theta.colSize();
    int numGhostCells = theta.ghostCellsPerBnd();
    Field rhsBuoyancyV(M, N + 1, 0.0, numGhostCells);

    for (int ii = theta.startIdxX(); ii <= theta.endIdxX(); ++ii)
    {
        for (int jj = theta.startIdxY(); jj <= theta.endIdxY(); ++jj)
        {
            rhsBuoyancyV(ii, jj) = raTimesPr*(theta(ii, jj) + theta(ii, jj-1))/2.0;
        }
    }

    return (rhsBuoyancyV);

}

#endif // COMPUTERHS_H
