#ifndef EXACTU
#define EXACTU

#include "Function2D.hpp"

/* This class specifies the exact solution u(x) in the PDE problem */

class ExactU: public Function2D {

  public:

    // Constructor
    ExactU();

    // Destructor
    ~ExactU();

    // u(x,T) 
    double evaluate(double x, double T);

};
#endif
