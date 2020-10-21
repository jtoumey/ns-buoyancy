/*----------------------------------------------------------------------------*\
  ME5311, Spring 2020
  Term Project
  Julian M. Toumey
--------------------------------------------------------------------------------
File
    field.hpp

Description
    A 2D array container for objects of type double. Inherits from Array2D, 
    leveraging constructors and member functions from that class. Also, adds 
    functionality for fields such as initialization and boundary conditions 
    specification.
\*----------------------------------------------------------------------------*/

#ifndef FIELD_H
#define FIELD_H

#include "array2d.hpp"

#include <fstream>
#include <time.h>
#include <cmath>

class Field : public Array2d  
{

protected:
    // Grid* activeGrid;

public:

    /// Inherit all constructors
    using Array2d::Array2d;

    /// Initialize the field with random numbers between 0 and 1
    void initializeRandom()
    {
        // Seed the RNG
        //srand(time(0));
        srand(0);

        // Fill the velocity field (physical cells) with random values
        for (int ii = startIdxX(); ii <= endIdxX(); ++ii)
        {
            for (int jj = startIdxY(); jj <= endIdxY(); ++jj)
            {
                // Create a random double between 0 and 1
                double randomNumx = (double)rand()/RAND_MAX; 
                arrayElements[ii][jj] = randomNumx;
            }
        }
    }

    /// DEBUG: Initialize the field with a uniform value 
    void initializeUniform(double uniformValue_)
    {
        // Set physical cells only
        for (int ii = startIdxX(); ii <= endIdxX(); ++ii)
        {
            for (int jj = startIdxY(); jj <= endIdxY(); ++jj)
            {
                arrayElements[ii][jj] = uniformValue_;
            }
        }        
    }

    /// DEBUG: Initialize the field with a trigonometric signal
    void initializeTrig(double Lx_, double Ly_)
    {
        // Because the field may have ghost cells, we need to decrement the 
        // running index ii by the numGhost so that the argument to the sin() 
        // start at the right value--zero
        int numGhost = (totalCellsx - numRows)/2;

        for (int ii = startIdxX(); ii <= endIdxX(); ++ii)
        {
            for (int jj = startIdxY(); jj <= endIdxY(); ++jj)
            {
                const double dx = Lx_/numRows;
                const double dy = Ly_/numCols;

                const double xx = dx*(ii - numGhost);
                const double yy = dy*(jj - numGhost);

                const double trigResult = std::cos(xx*(2*M_PI/Lx_))*std::cos(yy*(4*M_PI/Ly_));
                arrayElements[ii][jj] = trigResult;
            }
        }  
    }

    /// Set the West and East walls to a fixed value. Specific for u-velocity
    void setUBoundaryCondition(double boundaryValue_)
    {
        // Set boundary conditions
        for (int jj = startIdxY(); jj <= endIdxY(); ++jj)
        {
            // West wall
            arrayElements[startIdxX()][jj] = boundaryValue_;
            // East wall
            arrayElements[endIdxX()][jj] = boundaryValue_;
        }

        // Set South and North walls: No-slip
        int jjSouth = 0;
        int jjNorth = totalCellsy - 1;
        for (int ii = startIdxX(); ii <= endIdxX(); ++ii)
        {
            arrayElements[ii][jjSouth] = -arrayElements[ii][jjSouth+1];
            arrayElements[ii][jjNorth] = -arrayElements[ii][jjNorth-1];
        }
    }

    /// Set the North and South walls to a fixed value. Specific for v-velocity
    void setVBoundaryCondition(double boundaryValue_)
    {
        for (int ii = startIdxX(); ii <= endIdxX(); ++ii)
        {
            // South wall
            arrayElements[ii][startIdxY()] = boundaryValue_;
            // North wall
            arrayElements[ii][endIdxY()] = boundaryValue_;
        }

        // Set West and East walls: No-slip
        int iiWest = 0;
        int iiEast = totalCellsx - 1;
        for (int jj = startIdxY(); jj <= endIdxY(); ++jj)
        {
            arrayElements[iiWest][jj] = -arrayElements[iiWest+1][jj];
            arrayElements[iiEast][jj] = -arrayElements[iiEast-1][jj];
        }
    }

    /// Set the North and South walls to zero normal gradient by copying the 
    /// first off the wall cell values to the ghost cells. Specific to Theta in
    /// this case
    void setThetaBoundaryCondition()
    {
        // South wall: Set ghost cells to equal first-off-the-wall to enforce zero-
        // normal gradient BC
        int jjSouth = 0;
        // North wall: Same procedure
        int jjNorth = totalCellsy-1;

        for (int ii = startIdxX(); ii <= endIdxX(); ++ii)
        {
            arrayElements[ii][jjSouth] = arrayElements[ii][jjSouth+1];
            arrayElements[ii][jjNorth] = arrayElements[ii][jjNorth-1];
        }

        // West wall: 
        int iiWest = 0;
        // East wall:
        int iiEast = totalCellsx - 1;

        for (int jj = startIdxY(); jj <= endIdxY(); ++jj)
        {
            arrayElements[iiWest][jj] = 2 - arrayElements[iiWest+1][jj];
            arrayElements[iiEast][jj] = -arrayElements[iiEast-1][jj];
        }
    }

    /// Overload the + operator to allow field addition
    Field operator+(const Field& addend_)
    {
        // Construct a sum array with the same details as the addend
        Field sum
        (
            addend_.numRows, 
            addend_.numCols, 
            0.0, 
            addend_.numGhostCells
        );

        for (int ii = addend_.startIdxX(); ii <= addend_.endIdxX(); ++ii)
        {
            for (int jj = addend_.startIdxY(); jj <= addend_.endIdxY(); ++jj)
            {
                sum(ii, jj) = arrayElements[ii][jj] + addend_(ii, jj);
            }
        }
        return (sum);
    }

    /// Overload the - operator to allow field subtraction
    Field operator-(const Field& subtrahend_)
    {
        // Construct a difference array with the same details as the subtrahend
        Field difference
        (
            subtrahend_.numRows, 
            subtrahend_.numCols, 
            0.0, 
            subtrahend_.numGhostCells
        );

        for (int ii = subtrahend_.startIdxX(); ii <= subtrahend_.endIdxX(); ++ii)
        {
            for (int jj = subtrahend_.startIdxY(); jj <= subtrahend_.endIdxY(); ++jj)
            {
                difference(ii, jj) = arrayElements[ii][jj] - subtrahend_(ii, jj);
            }
        }
        return (difference);
    }

    /// Copy constructor
    Field(const Field& sourceArray_)
    {
        numRows = sourceArray_.rowSize();
        numCols = sourceArray_.colSize(); 

        totalCellsx = sourceArray_.totalCellsx;
        totalCellsy = sourceArray_.totalCellsy;

        physBndx = sourceArray_.physBndx;
        physBndy = sourceArray_.physBndy;

        numGhostCells = sourceArray_.numGhostCells;
        // Allocate memory on the heap
        arrayElements = new double*[totalCellsx];
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            arrayElements[ii] = new double[totalCellsy];
        }
        // Set each element to initial value
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            for (int jj = 0; jj < totalCellsy; ++jj)
            {
                arrayElements[ii][jj] = 0.0;
            }
        }
        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                arrayElements[ii][jj] = sourceArray_(ii, jj);
            }
        }
    }

    /// Overload the assignment operator to copy the internal array
    Field& operator=(const Field& sourceArray_)
    {
        // NOTE: This routine was added to avoid a memory leak. New allocation
        // was occurring in this routine without actually deleting previously
        // allocated arrays
        for (int ii = 0; ii < totalCellsx; ++ ii)
        {
            delete [] arrayElements[ii];
        }
        delete [] arrayElements;

        numRows = sourceArray_.rowSize();
        numCols = sourceArray_.colSize(); 

        totalCellsx = sourceArray_.totalCellsx;
        totalCellsy = sourceArray_.totalCellsy;

        physBndx = sourceArray_.physBndx;
        physBndy = sourceArray_.physBndy;

        numGhostCells = sourceArray_.numGhostCells;
        // Allocate memory on the heap
        arrayElements = new double*[totalCellsx];
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            arrayElements[ii] = new double[totalCellsy];
        }
        /////// FIXED 
        // Set each element to initial value
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            for (int jj = 0; jj < totalCellsy; ++jj)
            {
                arrayElements[ii][jj] = 0.0;
            }
        }

        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                arrayElements[ii][jj] = sourceArray_(ii, jj);
            }
        }
        return (*this);
    }

    /// Overload the assignment operator to copy the internal array from an 
    /// Array2D object
    void operator=(const Array2d& sourceArray_) 
    {
        numRows = sourceArray_.rowSize();
        numCols = sourceArray_.colSize(); 

        totalCellsx = sourceArray_.totalCellsX();
        totalCellsy = sourceArray_.totalCellsY();

        physBndx = sourceArray_.physBndX();
        physBndy = sourceArray_.physBndY();

        numGhostCells = sourceArray_.ghostCellsPerBnd();
        // Allocate memory on the heap
        arrayElements = new double*[totalCellsx];
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            arrayElements[ii] = new double[totalCellsy];
        }
        // Set each element to initial value
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            for (int jj = 0; jj < totalCellsy; ++jj)
            {
                arrayElements[ii][jj] = 0.0;
            }
        }
        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                arrayElements[ii][jj] = sourceArray_(ii, jj);
            }
        }
    }

    /// Copy constructor from an Array2d. Careful. Double careful.
    Field(const Array2d& sourceArray_)
    {
        numRows = sourceArray_.rowSize();
        numCols = sourceArray_.colSize(); 

        totalCellsx = sourceArray_.totalCellsX();
        totalCellsy = sourceArray_.totalCellsY();

        physBndx = sourceArray_.physBndX();
        physBndy = sourceArray_.physBndY();

        numGhostCells = sourceArray_.ghostCellsPerBnd();

        // Allocate memory on the heap
        arrayElements = new double*[totalCellsx];
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            arrayElements[ii] = new double[totalCellsy];
        }
        // Set each element to initial value
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            for (int jj = 0; jj < totalCellsy; ++jj)
            {
                arrayElements[ii][jj] = 0.0;
            }
        }
        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                arrayElements[ii][jj] = sourceArray_(ii, jj);
            }
        }
    }

    /// Overload the * operator for element-wise multiplication
    Field operator*(double multiplicand_) const
    {
        // Construct the product array with same details as the current array
        Field product
        (
            numRows, 
            numCols, 
            0.0,
            numGhostCells
        );

        for (int ii = startIdxX(); ii <= endIdxX(); ++ii)
        {
            for (int jj = startIdxY(); jj <= endIdxY(); ++jj)
            {
                product(ii, jj) = arrayElements[ii][jj]*multiplicand_;
            }
        }
        return (product);
    }
   
    void writeField(std::string fieldFileName)
    {
        std::ofstream fieldOut;
        fieldOut.open(fieldFileName);

        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                fieldOut << std::setw(20) << std::left << arrayElements[ii][jj];
            }
            fieldOut << std::endl;
        }
    }
};

#endif // FIELD_H