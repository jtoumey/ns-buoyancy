/*----------------------------------------------------------------------------*\
  ME5311, Spring 2020
  Term Project
  Julian M. Toumey
--------------------------------------------------------------------------------
File
    advanceRungeKutta.hpp

Description
    Based on the state vector, the derivative of the state vector, and a `step`
    flag, this function advances the state per the Runge-Kutta routine of 
    Spalart, Moses, and Rogers.
\*----------------------------------------------------------------------------*/

#ifndef ADVANCERUNGEKUTTA_H
#define ADVANCERUNGEKUTTA_H

#include "field.hpp"

Field advanceRungeKutta(Field& u0, Field& rhs0, Field& rhs1, int step, double dt)
{
    double rkCoeff1, rkCoeff2;

    switch (step) {
        case 1 :
            rkCoeff1 = 8.0/15.0; 
            rkCoeff2 = 0.0;
            break;
        case 2 :
            rkCoeff1 = -17.0/60.0; 
            rkCoeff2 = 5.0/12.0;
            break;
        case 3 :
            rkCoeff1 = -5.0/12.0;
            rkCoeff2 = 3.0/4.0;
            break;
        default :
            std::cout << "Invalid RK step." << std::endl;
            exit(0);
    }

    Field advancedField = u0 + rhs0*rkCoeff1*dt + rhs1*rkCoeff2*dt;

    return advancedField;
}

#endif // ADVANCERUNGEKUTTA_H