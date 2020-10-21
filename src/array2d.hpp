/*----------------------------------------------------------------------------*\
  ME5311, Spring 2020
  Term Project
  Julian M. Toumey
--------------------------------------------------------------------------------
File
    array2d.hpp
    
Description
    A 2D array container for objects of type double. 
    Storage is allocated on the heap when anything other than the default 
    constructor is called.
\*----------------------------------------------------------------------------*/

#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

/// Support struct for holding array bounds, either physical or logical
struct BoundIndex
{
    int startIndex;
    int endIndex;

    /// Default constructor
    BoundIndex()
    {
        startIndex = 0;
        endIndex = 0;
    }

    /// Constructor with given bounds
    BoundIndex(int startIndex_, int endIndex_)
    {
        startIndex = startIndex_;
        endIndex = endIndex_;
    }
};

class Array2d
{
protected:

    /// Physical rows and columns in the array per coordinate direction
    int numRows;
    int numCols;

    /// Total cells in the array (including ghost cells) per coordinate 
    /// direction
    int totalCellsx;
    int totalCellsy;

    // Ghost cells PER BOUNDARY in a given coordinate direction
    int numGhostCells;

    /// Physical and logical boundary structures
    BoundIndex physBndx;
    BoundIndex physBndy;
    BoundIndex logicalBndx;
    BoundIndex logicalBndy;

    /// Actual array elements: Pointer to an array of pointers to heap storage
    double** arrayElements;

public:

    // * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //
    
    /// Default constructor
    Array2d()
    {
        numRows = 0;
        numCols = 0;
    }

    /// Constructor with physical dimensions only
    Array2d(int numRows_, int numCols_)
    {
        numRows = numRows_;
        numCols = numCols_;

        totalCellsx = numRows_;
        totalCellsy = numCols_;

        numGhostCells = 0;

        physBndx.startIndex = 0;
        physBndx.endIndex = numRows_ - 1;

        physBndy.startIndex = 0;
        physBndy.endIndex = numCols_ - 1;

        // Allocate memory on heap
        arrayElements = new double*[totalCellsx];
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            arrayElements[ii] = new double[totalCellsy];
        }
    }

    /// Constructor with physical dimensions and initial value
    Array2d(int numRows_, int numCols_, double initialValue_)
    {
        numRows = numRows_;
        numCols = numCols_;

        totalCellsx = numRows_;
        totalCellsy = numCols_;

        numGhostCells = 0;

        physBndx.startIndex = 0;
        physBndx.endIndex = numRows_ - 1;

        physBndy.startIndex = 0;
        physBndy.endIndex = numCols_ - 1;

        // Allocate memory on the heap
        arrayElements = new double*[totalCellsx];
        for (int ii = 0; ii < numRows; ++ii)
        {
            arrayElements[ii] = new double[totalCellsy];
        }

        // Set each element to zero
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            for (int jj = 0; jj < totalCellsy; ++jj)
            {
                arrayElements[ii][jj] = 0.0;
            }
        }
    }

    /// Constructor with physical dimensions, initial value, and number of ghost 
    /// cells
    Array2d(int numRows_, int numCols_, double initialValue_, 
            int numGhostCells_)
    {
        numRows = numRows_;
        numCols = numCols_;

        physBndx.startIndex = numGhostCells_;
        physBndx.endIndex = numRows_ + numGhostCells_ - 1;

        physBndy.startIndex = numGhostCells_;
        physBndy.endIndex = numCols_ + numGhostCells_ - 1;

        numGhostCells = numGhostCells_;
     
        totalCellsx = numRows_ + 2*numGhostCells_;
        totalCellsy = numCols_ + 2*numGhostCells_;

        // Allocate memory on the heap
        arrayElements = new double*[totalCellsx];
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            arrayElements[ii] = new double[totalCellsy];
        }

        // Set each element to 0.0
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            for (int jj = 0; jj < totalCellsy; ++jj)
            {
                arrayElements[ii][jj] = 0.0;
            }
        }
        // Set interior elements to initial value
        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                arrayElements[ii][jj] = initialValue_;
            }
        }
    }


    // * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * * //

    /// Free the dynamically allocated array
    virtual ~Array2d()
    {
        for (int ii = 0; ii < totalCellsx; ++ii) 
        {
            delete [] arrayElements[ii];
        }
        delete [] arrayElements;
    }


    // * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

    /// Verify the array access bounds are valid 
    void checkBounds(int rowIdx_, int colIdx_) const
    {
        // Note that we are just checking that indices are within logical bounds

        if (rowIdx_ < 0 || rowIdx_ > totalCellsx - 1)
        {
            std::cout << "ERROR: Row index " << rowIdx_ 
                      << " is out of allowed bounds [" << 0 << ", " 
                      << totalCellsx - 1 << "]" << std::endl;
            exit(0);
        }
        if (colIdx_ < 0 || colIdx_ > totalCellsy - 1)
        {
            std::cout << "ERROR: Col index " << colIdx_ 
                      << " is out of allowed bounds [" << 0 << ", " 
                      << totalCellsy - 1 << "]" << std::endl;
            exit(0);
        }
    }

    double max()
    {
        double maxValue = 0.0;
        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                if (arrayElements[ii][jj] > maxValue)
                {
                    maxValue = arrayElements[ii][jj];
                }
            }
        }
        return (maxValue);
    }

    /// Return the maximum element in the array
    double maxAbs()
    {
        double maxAbsValue = 0.0;

        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                if (std::abs(arrayElements[ii][jj]) > maxAbsValue)
                {
                    maxAbsValue = std::abs(arrayElements[ii][jj]);
                }
            }
        }
        return (maxAbsValue);
    }

    /// Return the array with each term squared
    Array2d square()
    {
        // Construct a squared array with the same details as the current array
        Array2d arraySquared
        (
            numRows, 
            numCols, 
            0.0, 
            numGhostCells
        );

        // Square each term and save to the output array
        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                arraySquared(ii, jj) = arrayElements[ii][jj]*arrayElements[ii][jj];
            }
        }
        return (arraySquared);
    }

    /// Return the sum of all terms in the array
    double sum()
    {
        double arraySum = 0.0;

        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                arraySum += arrayElements[ii][jj];
            }
        }
        return (arraySum);
    }

    /// Return the array row size (physical size)
    int rowSize() const
    {
        return (numRows);
    }

    /// Return the array column size (physical size)
    int colSize() const
    {
        return (numCols);
    }

    /// Return the physical start index for coordinate direction x
    int startIdxX() const
    {
        return (physBndx.startIndex);
    }

    /// Return the physical start index for coordinate direction y
    int startIdxY() const
    {
        return (physBndy.startIndex);
    }

    /// Return the physical end index for coordinate direction x
    int endIdxX() const
    {
        return (physBndx.endIndex);
    }

    /// Return the physical end index for coordinate direction y
    int endIdxY() const
    {
        return (physBndy.endIndex);
    }

    /// Return the total number of cells (physical and logical) in x
    int totalCellsX() const
    {
        return (totalCellsx);
    }

    /// Return the total number of cells (physical and logical) in y
    int totalCellsY() const
    {
        return (totalCellsy);
    }

    /// Return the number of ghost cells per boundary (constant per boundary)
    int ghostCellsPerBnd() const
    {
        return (numGhostCells);
    }

    /// Return the bound indices for coordinate direction x
    BoundIndex physBndX() const
    {
        return (physBndx);
    }

    /// Return the bound indices for coordinate direction y
    BoundIndex physBndY() const
    {
        return (physBndy);
    }


    // * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //

    /// Overload the () operator for index-wise element access
    double& operator()(int rowIdx, int colIdx) const
    {
        // Verify the supplied indices
        checkBounds(rowIdx, colIdx);

        // Dereference the array of elements
        return (arrayElements[rowIdx][colIdx]);
    }

    /// 
    Array2d& operator+=(const Array2d& addend_)
    {
        // Array2d sum
        // (
        //     addend_.numRows, 
        //     addend_.numCols, 
        //     0.0, 
        //     addend_.numGhostCells
        // );
        // Loop and add term by term
        for (int ii = addend_.startIdxX(); ii <= addend_.endIdxX(); ++ii)
        {
            for (int jj = addend_.startIdxY(); jj <= addend_.endIdxY(); ++jj)
            {
                // sum(ii, jj) = arrayElements[ii][jj] + addend_(ii, jj);
                arrayElements[ii][jj] = arrayElements[ii][jj] + addend_(ii, jj);
            }
        }
        return (*this);        
    }

    /// Overload the + operator to allow field addition
    Array2d& operator+(const Array2d& addend_)
    {
        Array2d sum(*this);
        return sum += addend_;
        // Construct a sum array with the same details as the addend
        // Array2d sum
        // (
        //     addend_.numRows, 
        //     addend_.numCols, 
        //     0.0, 
        //     addend_.numGhostCells
        // );

        // // Loop and add term by term
        // for (int ii = addend_.startIdxX(); ii <= addend_.endIdxX(); ++ii)
        // {
        //     for (int jj = addend_.startIdxY(); jj <= addend_.endIdxY(); ++jj)
        //     {
        //         sum(ii, jj) = arrayElements[ii][jj] + addend_(ii, jj);
        //     }
        // }
        // return (sum);

    }

    /// Overload the - operator to allow field subtraction
    Array2d operator-(const Array2d& subtrahend_)
    {
        // Construct a difference array with the same details as the subtrahend
        Array2d difference
        (
            subtrahend_.numRows, 
            subtrahend_.numCols, 
            0.0, 
            subtrahend_.numGhostCells
        );

        // Loop and subtract term bu term
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
    Array2d(const Array2d& sourceArray_)
    {
        numRows = sourceArray_.numRows;
        numCols = sourceArray_.numCols;

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

        // Set all elements to zero
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            for (int jj = 0; jj < totalCellsy; ++jj)
            {
                arrayElements[ii][jj] = 0.0;
            }
        }

        // Copy over physical elements
        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                arrayElements[ii][jj] = sourceArray_(ii, jj);
            }
        }
    }

    /// Overload the assignment operator to copy the internal array
    Array2d& operator=(const Array2d& sourceArray_)
    {
        // NOTE: This routine was added to avoid a memory leak. New allocation
        // was occurring in this routine without actually deleting previously
        // allocated arrays
        for (int ii = 0; ii < totalCellsx; ++ ii)
        {
            delete [] arrayElements[ii];
        }
        delete [] arrayElements;
        
        numRows = sourceArray_.numRows;
        numCols = sourceArray_.numCols;
        
        totalCellsx = sourceArray_.totalCellsx;
        totalCellsy = sourceArray_.totalCellsy;

        physBndx = sourceArray_.physBndx;
        physBndy = sourceArray_.physBndy;

        // Allocate memory on the heap
        arrayElements = new double*[totalCellsx];
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            arrayElements[ii] = new double[totalCellsy];
        }

        // Set all elements to zero
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            for (int jj = 0; jj < totalCellsy; ++jj)
            {
                arrayElements[ii][jj] = 0.0;
            }
        }

        // Copy over physical elements
        for (int ii = physBndx.startIndex; ii <= physBndx.endIndex; ++ii)
        {
            for (int jj = physBndy.startIndex; jj <= physBndy.endIndex; ++jj)
            {
                arrayElements[ii][jj] = sourceArray_(ii, jj);
            }
        }
        return (*this);
    }

    /// Overload the * operator for element-wise multiplication
    Array2d operator*(double multiplicand_) const
    {
        // Construct the product array with same details as the current array
        Array2d product
        (
            numRows, 
            numCols, 
            0.0,
            numGhostCells
        );

        // Loop multiplying term by term
        for (int ii = startIdxX(); ii <= endIdxX(); ++ii)
        {
            for (int jj = startIdxY(); jj <= endIdxY(); ++jj)
            {
                product(ii, jj) = arrayElements[ii][jj]*multiplicand_;
            }
        }
        return (product);
    }

    /// DEBUG: Print array contents with string label
    void printArray(std::string arrayLabel)
    {
        std::cout << "====  * " << arrayLabel << " *  ====" << std::endl;
        const int lineWidth = 8;
        std::cout << std::setprecision(4) << std::fixed;
        for (int ii = 0; ii < totalCellsx; ++ii)
        {
            for (int jj = 0; jj < totalCellsy; ++jj)
            {
                std::cout << std::setw(lineWidth) << std::right 
                          << arrayElements[ii][jj] << "\t";
            }
            std::cout << std::endl;
        }
    }

    /// DEBUG: Print array dimensions with string label
    void printBounds(std::string arrayLabel) const
    {
        std::cout << "====  + " << arrayLabel << " +  ====" << std::endl;
        std::cout << "Cell count: " << std::endl;
        std::cout << "  [" << totalCellsx << ", " << totalCellsy << "]" 
                  << std::endl;
        std::cout << "Physical array bounds: " << std::endl;
        std::cout << "  x: [" << physBndx.startIndex << ", " 
                  << physBndx.endIndex << "]" << std::endl;
        std::cout << "  y: [" << physBndy.startIndex << ", " 
                  << physBndy.endIndex << "]" << std::endl;
    }

};

#endif // ARRAY2D_H