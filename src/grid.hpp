/*----------------------------------------------------------------------------*\
  ME5311, Spring 2020
  Term Project
  Julian M. Toumey
--------------------------------------------------------------------------------
File
    grid.hpp

Description
    A container for the grid details, including node locations.

\*----------------------------------------------------------------------------*/

#ifndef GRID_H
#define GRID_H

#include <fstream>
#include <iomanip>
#include <string>


class Grid
{
private:

    /// Grid details 
    int numCellsx;
    int numCellsy;

    double lengthx;
    double lengthy;

    double dx;
    double dy;

    /// Arrays of node location details
    double* cellCentersx;
    double* cellCentersy;

    double* nodeCoordinatesx;
    double* nodeCoordinatesy;

public:

    // * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //
    
    /// Constructor with number of cells and size
    Grid(int numCellsx_, int numCellsy_, double lengthx_, double lengthy_) :
         numCellsx(numCellsx_), numCellsy(numCellsy_), lengthx(lengthx_),
         lengthy(lengthy_)
    {
        double dx = (double) lengthx/numCellsx;
        double dy = (double) lengthy/numCellsy;

        // Form cell center coordinates
        cellCentersx = new double[numCellsx];

        for (int ii = 0; ii < numCellsx; ++ii)
        {
            double xc = dx/2.0 + ii*dx;
            cellCentersx[ii] = xc;
        }

        cellCentersy = new double[numCellsy];

        for (int ii = 0; ii < numCellsy; ++ii)
        {
            double yc = dy/2.0 + ii*dy;
            cellCentersy[ii] = yc;
        }

        // Form node coordinates
        nodeCoordinatesx = new double[numCellsx + 1];
        
        for (int ii = 0; ii < numCellsx + 1; ++ii)
        {
            nodeCoordinatesx[ii] = ii*dx;
        }

        nodeCoordinatesy = new double[numCellsy + 1];
        
        for (int ii = 0; ii < numCellsy + 1; ++ii)
        {
            nodeCoordinatesy[ii] = ii*dy;
        }
    }

    // * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //
    ~Grid()
    {
        delete [] cellCentersx;
        delete [] cellCentersy;
        delete [] nodeCoordinatesx;
        delete [] nodeCoordinatesy;
    } 

    /// Write grid coordinates to output files
    void writeGrid()
    {
        // Write x node locations
        std::ofstream xNodeLoc;
        xNodeLoc.open("nodeLocation_x.dat");
        for (int ii = 0; ii < numCellsx + 1; ++ii)
        {
            for (int jj = 0; jj < numCellsy + 1; ++jj)
            {
                xNodeLoc << std::setw(20) << std::left << nodeCoordinatesx[ii];
            }
            xNodeLoc << std::endl;
        }
        // Write y node locations
        std::ofstream yNodeLoc;
        yNodeLoc.open("nodeLocation_y.dat");
        for (int ii = 0; ii < numCellsx + 1; ++ii)
        {
            for (int jj = 0; jj < numCellsy + 1; ++jj)
            {
                yNodeLoc << std::setw(20) << std::left << nodeCoordinatesy[jj];
            }
            yNodeLoc << std::endl;
        }        
    }
    
};

#endif // GRID_H