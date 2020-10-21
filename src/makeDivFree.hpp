

#ifndef MAKEDIVFREE_H
#define MAKEDIVFREE_H

void makeDivFree(Field& uVel, Field& vVel, double dx, double dy)
{
    // 
    int M = vVel.rowSize();
    int N = uVel.colSize();
    
    // Solve the Poisson equation to get div-free velocity
    Array2d divVelocity = div(uVel, vVel, dx, dy);
    Array2d pressure(M, N, 0.0);
    pressure = solvePoisson(divVelocity, dx, dy);

    // Compute pressure gradient: dp/dx
    Array2d dpDx = gradx(pressure, dx);
    // dp/dy
    Array2d dpDy = grady(pressure, dy);
    
    // Add pressure gradient to velocity
    uVel = uVel - dpDx;
    vVel = vVel - dpDy;

}

#endif // MAKEDIVFREE_H
