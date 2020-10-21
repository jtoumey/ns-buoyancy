/*----------------------------------------------------------------------------*\
  ME5311, Spring 2020
  Term Project
  Julian M. Toumey
--------------------------------------------------------------------------------
File
    fieldCalculus.hpp

Description
    Collection of functions for computing the discretized convection terms in 
    the Navier-Stokes equations. 
\*----------------------------------------------------------------------------*/

#ifndef FIELDCALCULUS_H
#define FIELDCALCULUS_H

#include "array2d.hpp"
#include "field.hpp"

Array2d div(Field& uVel, Field& vVel, double dx, double dy)
{
    // v-Velocity is stored at center of face in x-dir
    int M = vVel.rowSize();
    // u-Velocity is stored at center of face in y-dir
    int N = uVel.colSize();

    Array2d divVelocity(M, N, 0.0);
    for (int ii = uVel.startIdxX(); ii <= uVel.endIdxX() - 1; ++ii)
    {
        for (int jj = vVel.startIdxY(); jj <= vVel.endIdxY() - 1; ++jj)
        {
            // No ghost cells in our velocity divergence array, so we set the 
            // indices back by one
            divVelocity(ii - 1, jj - 1) = (uVel(ii+1, jj) - uVel(ii, jj))/dx + 
                                          (vVel(ii, jj+1) - vVel(ii, jj))/dy;
        }
    }
    return (divVelocity);
}

Array2d gradx(Array2d& inputField, double dx)
{
    int M = inputField.rowSize();
    int N = inputField.colSize();

    // TODO: Fix this ad-hoc statement
    int numGhostCells = 1;
    // Problem is: The input array has no ghost cells, but I want to pad out the
    // output array with one ghost cell. 
    // int numGhostCells = inputField.ghostCellsPerBnd();

    Array2d gradField(M + 1, N, 0.0, numGhostCells);
    for (int ii = inputField.startIdxX(); ii <= inputField.endIdxX() - 1; ++ii)
    {
        for (int jj = inputField.startIdxY(); jj <= inputField.endIdxY(); ++jj)
        {
            // Unfortunate indexing trick to get the pressure array (cell 
            // centered, no ghost cells) into an array directly compatible with
            // the uVelocity
            gradField(ii + 2, jj + 1) = (inputField(ii + 1, jj) - 
                                         inputField(ii, jj))/dx;
        }
    }
    return (gradField);
}

Array2d grady(Array2d& inputField, double dy)
{
    int M = inputField.rowSize();
    int N = inputField.colSize();

    // TODO: Fix this ad-hoc statement
    int numGhostCells = 1;

    Array2d gradField(M, N + 1, 0.0, numGhostCells);
    for (int ii = inputField.startIdxX(); ii <= inputField.endIdxX(); ++ii)
    {
        for (int jj = inputField.startIdxY(); jj <= inputField.endIdxY() - 1; ++jj)
        {
            gradField(ii + 1, jj + 2) = (inputField(ii, jj + 1) - 
                                         inputField(ii, jj))/dy;
        }
    }
    return (gradField);
}

Field laplacianx(Field& uVel, double dx, double dy)
{
    int M = uVel.rowSize();
    int N = uVel.colSize();
    int ngc = uVel.ghostCellsPerBnd();
    Field laplacianField(M, N, 0.0, ngc);

    for (int ii = uVel.startIdxX(); ii <= uVel.endIdxX(); ++ii)
    {
        for (int jj = uVel.startIdxY(); jj <= uVel.endIdxY(); ++jj)
        {   
            double lpfx = (uVel(ii-1, jj) - 2.0*uVel(ii, jj) + uVel(ii+1, jj))/(dx*dx);
            double lpfy = (uVel(ii, jj-1) - 2.0*uVel(ii, jj) + uVel(ii, jj+1))/(dy*dy);
            laplacianField(ii, jj) = lpfx + lpfy;
        }
    }
    return (laplacianField);
}

Field laplaciany(Field& vVel, double dx, double dy)
{
    int M = vVel.rowSize();
    int N = vVel.colSize();
    int ngc = vVel.ghostCellsPerBnd();
    Field laplacianField(M, N, 0.0, ngc);

    for (int ii = vVel.startIdxX(); ii <= vVel.endIdxX(); ++ii)
    {
        for (int jj = vVel.startIdxY(); jj <= vVel.endIdxY(); ++jj)
        {
            double lpfx = (vVel(ii-1, jj) - 2.0*vVel(ii, jj) + vVel(ii+1, jj))/(dx*dx);
            double lpfy = (vVel(ii, jj-1) - 2.0*vVel(ii, jj) + vVel(ii, jj+1))/(dy*dy);
            laplacianField(ii, jj) = lpfx + lpfy;
        }
    }
    return (laplacianField);
}

Field laplacianTheta(Field& theta, double dx, double dy)
{
    int M = theta.rowSize();
    int N = theta.colSize();
    int ngc = theta.ghostCellsPerBnd();
    Field laplacianField(M, N, 0.0, ngc);

    for (int ii = theta.startIdxX(); ii <= theta.endIdxX(); ++ii)
    {
        for (int jj = theta.startIdxY(); jj <= theta.endIdxY(); ++jj)
        {
            double lpfx = (theta(ii-1, jj) - 2.0*theta(ii, jj) + theta(ii+1, jj))/(dx*dx);
            double lpfy = (theta(ii, jj-1) - 2.0*theta(ii, jj) + theta(ii, jj+1))/(dy*dy);
            laplacianField(ii, jj) = lpfx + lpfy;
        }
    }
    return (laplacianField);    
}

#endif // FIELDCALCULUS_H