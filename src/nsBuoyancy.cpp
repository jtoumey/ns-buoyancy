/*----------------------------------------------------------------------------*\
  ME5311, Spring 2020
  Term Project
  Julian M. Toumey
--------------------------------------------------------------------------------
File
    nsBuoyancy.cpp

Description
    Main driver file for solving the non-dimensional Navier-Stokes equations on
    a 2D structured Cartesian grid generated ad hoc.
    Includes convective and diffusive fluxes, as well as buoyancy for vertical
    y-direction.

    Requires compiler support for C++17. 
    Compile via `g++ nsBuoyancy.cpp -std=c++17` or similar.
\*----------------------------------------------------------------------------*/


#include <cstdlib>

#include "grid.hpp"
#include "array2d.hpp"
#include "field.hpp"
#include "fieldCalculus.hpp"
#include "solvePoisson.hpp"
#include "computeRhs.hpp"
#include "makeDivFree.hpp"
#include "advanceRungeKutta.hpp"

int main (int argc, char** argv)
{
    /// These parameters are set to default values. If command line args are 
    /// supplied, they may be overwritten
    // Grid Setup 
    int M = 64;
    int N = 128;

    // Simulation parameters 
    const double prandtlNum = 1.0;
    double rayleighNum = 5.0e8;

    bool variableTimeStepFlag = false;
    // Now, see if any command line args were supplied and overwrite the default
    // setup
    if (argc == 1)
    {
        // Default configuration with no arguments. Proceed as normal.
    }
    else if (argc == 4)
    {
        // Parse numCellsx, numCellsy, and Rayleigh number
        M = atoi(argv[1]);
        N = atoi(argv[2]);
        rayleighNum = atof(argv[3]);
    }
    else if (argc == 5)
    {
        // Parse numCellsx, numCellsy, Rayleigh number, varDt
        M = atoi(argv[1]);
        N = atoi(argv[2]);
        rayleighNum = atof(argv[3]);
        variableTimeStepFlag = true;
    }
    else
    {
        std::cout << "ERROR: When supplying command line arguments, you must"
                  << " enter \n\n"
                  << "  ./a.out M N RayleighNum (varDt)\n" 
                  << "\n  M: Number of cells in x "
                  << "\n  N: Number of cells in y "
                  << "\n  RayleighNum: Rayleigh number" 
                  << "\n  (varDt): OPTIONAL. Enable variable time-stepping."
                  << std::endl;
        exit(0);
    }
    

    double Lx = 0.5;
    double Ly = 1.0;
    double dx = Lx/static_cast<double> (M);
    double dy = Ly/static_cast<double> (N);
    int numGhostCells = 1;

    Grid activeGrid(M, N, Lx, Ly);
    activeGrid.writeGrid();

    // Initial field setup
    Field uVel(M + 1, N, 0.0, numGhostCells);
    // uVel.initializeRandom(); 
    uVel.initializeUniform(0.0);
    uVel.setUBoundaryCondition(0.0);
   
    Field vVel(M, N + 1, 0.0, numGhostCells);
    // vVel.initializeRandom();
    vVel.initializeUniform(0.0);
    vVel.setVBoundaryCondition(0.0);

    Field theta(M, N, 0.0, numGhostCells);
    // theta.initializeRandom();
    theta.initializeUniform(0.5);
    // theta.initializeTrig(Lx, Ly);
    theta.setThetaBoundaryCondition();

    // Make velocity divergence-free
    makeDivFree(uVel, vVel, dx, dy);
    uVel.setUBoundaryCondition(0.0);
    vVel.setVBoundaryCondition(0.0);


    // Time integration parameters
    double cflTarget = 1.0;
    int numTimeSteps = 5000;
    double t0 = 0.0;
    double tEnd;
    double tCurrent = t0;
    double dt = 1e-6;

    // Save the max divergence for stability checking
    double maxDiv = 0.0;
    // t_{n - 1} divergence--for comparison
    double maxDiv0 = 0.0;
    double maxVel = 0.0;
    double maxVel0 = 0.0;

   
    // Write initial fields to disk
    std::cout << std::setprecision(8) << std::fixed;   
    std::string tCurr = std::to_string(tCurrent);
    theta.writeField("theta_t" + tCurr + ".dat");
    uVel.writeField("uVel_t" + tCurr + ".dat");
    vVel.writeField("vVel_t" + tCurr + ".dat");

    int tt = 0;
    std::cout << "Starting time integration..." << std::endl;
    // for (int tt = 0; tt < numTimeSteps; ++tt)
    while (tCurrent < 1e-3)
    {
        // Update current physical time
        tCurrent = t0 + (tt + 1)*dt; 

        // * * * * Runge-Kutta loop * * * * //
        // * * RK sub-step 1 * * //
        // Theta
        Field rhsViscTheta = computeRhsViscTheta(theta, dx, dy, prandtlNum);
        Field rhsConvectTheta = computeRhsTheta(theta, uVel, vVel, dx, dy);
        Field rhsTheta = rhsConvectTheta + rhsViscTheta;
        Field thetaAlpha = advanceRungeKutta(theta, rhsTheta, rhsTheta, 1, dt);
        thetaAlpha.setThetaBoundaryCondition();

        // u-velocity
        Field rhsConvectU = computeRhsU(uVel, vVel, dx, dy);
        Field rhsViscU = computeRhsViscU(uVel, dx, dy, prandtlNum);
        Field rhsU = rhsConvectU + rhsViscU;
        Field uVelAlpha = advanceRungeKutta(uVel, rhsU, rhsU, 1, dt);
        uVelAlpha.setUBoundaryCondition(0.0);

        // v-velocity
        Field rhsConvectV = computeRhsV(uVel, vVel, dx, dy);
        Field rhsViscV = computeRhsViscV(vVel, dx, dy, prandtlNum);
        Field rhsBuoyancyV = computeRhsBuoyancyV(theta, rayleighNum, prandtlNum);
        Field rhsV = rhsConvectV + rhsBuoyancyV + rhsViscV;
        Field vVelAlpha = advanceRungeKutta(vVel, rhsV, rhsV, 1, dt);
        vVelAlpha.setVBoundaryCondition(0.0);

        makeDivFree(uVelAlpha, vVelAlpha, dx, dy);
        uVelAlpha.setUBoundaryCondition(0.0);
        vVelAlpha.setVBoundaryCondition(0.0);


        // * * RK sub-step 2 * * //
        Field rhsViscThetaAlpha = computeRhsViscTheta(thetaAlpha, dx, dy, prandtlNum);
        Field rhsConvectThetaAlpha = computeRhsTheta(thetaAlpha, uVelAlpha, vVelAlpha, dx, dy);
        Field rhsThetaAlpha = rhsConvectThetaAlpha + rhsViscThetaAlpha;
        Field thetaAlpha1 = advanceRungeKutta(thetaAlpha, rhsTheta, rhsThetaAlpha, 2, dt);

        Field rhsConvectUalpha = computeRhsU(uVelAlpha, vVelAlpha, dx, dy);
        Field rhsViscUalpha = computeRhsViscU(uVelAlpha, dx, dy, prandtlNum);
        Field rhsUalpha = rhsConvectUalpha + rhsViscUalpha;
        Field uVelAlpha1 = advanceRungeKutta(uVelAlpha, rhsU, rhsUalpha, 2, dt);

        Field rhsConvectValpha = computeRhsV(uVelAlpha, vVelAlpha, dx, dy);
        Field rhsViscValpha = computeRhsViscV(vVelAlpha, dx, dy, prandtlNum);
        Field rhsBuoyancyValpha = computeRhsBuoyancyV(thetaAlpha, rayleighNum, prandtlNum);
        Field rhsValpha = rhsConvectValpha + rhsViscValpha + rhsBuoyancyValpha;
        Field vVelAlpha1 = advanceRungeKutta(vVelAlpha, rhsV, rhsValpha, 2, dt);

        thetaAlpha1.setThetaBoundaryCondition();
        uVelAlpha1.setUBoundaryCondition(0.0);
        vVelAlpha1.setVBoundaryCondition(0.0);
        makeDivFree(uVelAlpha1, vVelAlpha1, dx, dy);
        uVelAlpha1.setUBoundaryCondition(0.0);
        vVelAlpha1.setVBoundaryCondition(0.0);

        // * * RK sub-step 3 * * //
        Field rhsViscThetaAlpha1 = computeRhsViscTheta(thetaAlpha1, dx, dy, prandtlNum);
        Field rhsConvectThetaAlpha1 = computeRhsTheta(thetaAlpha1, uVelAlpha1, vVelAlpha1, dx, dy);
        Field rhsThetaAlpha1 = rhsConvectThetaAlpha1 + rhsViscThetaAlpha1;
        Field thetaNp1 = advanceRungeKutta(thetaAlpha1, rhsThetaAlpha, rhsThetaAlpha1, 3, dt);
        theta = thetaNp1;

        Field rhsConvectUalpha1 = computeRhsU(uVelAlpha1, vVelAlpha1, dx, dy);
        Field rhsViscUalpha1 = computeRhsViscU(uVelAlpha1, dx, dy, prandtlNum);
        Field rhsUalpha1 = rhsConvectUalpha1 + rhsViscUalpha1;
        Field uVelNp1 = advanceRungeKutta(uVelAlpha1, rhsUalpha, rhsUalpha1, 3, dt);
        
        Field rhsConvectValpha1 = computeRhsV(uVelAlpha1, vVelAlpha1, dx, dy);
        Field rhsViscValpha1 = computeRhsViscV(vVelAlpha1, dx, dy, prandtlNum);
        Field rhsBuoyancyValpha1 = computeRhsBuoyancyV(thetaAlpha1, rayleighNum, prandtlNum);
        Field rhsValpha1 = rhsConvectValpha1 + rhsViscValpha1 + rhsBuoyancyValpha1;
        Field vVelNp1 = advanceRungeKutta(vVelAlpha1, rhsValpha, rhsValpha1, 3, dt);

        theta.setThetaBoundaryCondition();
        uVelNp1.setUBoundaryCondition(0.0);
        vVelNp1.setVBoundaryCondition(0.0);
        // uVelNp1.printArray("uVel Np1::blowing up");
        makeDivFree(uVelNp1, vVelNp1, dx, dy);
        uVelNp1.setUBoundaryCondition(0.0);
        vVelNp1.setVBoundaryCondition(0.0);

        uVel = uVelNp1;
        vVel = vVelNp1;
        uVel.setUBoundaryCondition(0.0);
        vVel.setVBoundaryCondition(0.0);

    
        // Compute maximum CFL and new dt based on this
        double maxCflc = uVel.maxAbs()*dt/dx + vVel.maxAbs()*dt/dx;
        double maxCfld = prandtlNum*dt/(dx*dx);
        double dtNew;
        if (variableTimeStepFlag)
        {
            double maxCfl = 0;
            if (maxCflc > maxCfld)
            {
                maxCfl = maxCflc;
            }
            else
            {
                maxCfl = maxCfld;
            }
            dtNew = (0.5*cflTarget/maxCfl + 0.5)*dt*0.1;
        }
        else
        {
            dtNew = dt;
        }

        maxVel0 = maxVel;
        maxVel = uVel.max();
        maxDiv = (div(uVel, vVel, dx, dy)).maxAbs();
        
        if (tt > 5)
        {
            double diffVel = std::abs(maxVel - maxVel0)/maxVel0;
            // std:: cout << "\nChange in max uVel: " << diffVel << std::endl;
            if (diffVel > 0.1)
            // Test for instability
            {
                std::cout << "\nWARNING: Velocity instability detected.\n" 
                          << "\tCutting time-step by factor of 10.\n"
                          << std::endl;
                dtNew = dt*0.1;
            }
        }

        // Screen output: Time-step
        std::cout << "\nTime-step: " << tt + 1 << "\t\tCurrent time: " 
                  << tCurrent << "\t Max divergence: " << std::scientific
                  << maxDiv << std::endl;

        std::cout << "dt: " << dt << "\tMax CFL_c: " << maxCflc 
                  << "\t Max CFL_d: " << maxCfld << "\tNew dt: " 
                  << dtNew  << std::endl;

        tEnd = tCurrent;
        
        // After screen output, update the time-step size
        dt = dtNew;

        // Dump output 
        if (tt % 10 == 0)
        {
            std::string tCurr = std::to_string(tCurrent);

            theta.writeField("theta_t" + tCurr + ".dat");
            uVel.writeField("uVel_t" + tCurr + ".dat");
            vVel.writeField("vVel_t" + tCurr + ".dat");
        }

        // Detect instability and exit clean
        if ((div(uVel, vVel, dx, dy)).maxAbs() > 1e8)
        {
            std::cout << "ERROR: Velocity instability detected." << std::endl;
            exit (0);
        }
        ++tt;
    } 
}
