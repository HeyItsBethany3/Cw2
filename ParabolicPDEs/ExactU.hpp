#ifndef EXACTU
#define EXACTU

#include "Function2D.hpp"

// This function specifies the exact solution for u in our PDE problem
class ExactU: public Function2D {

  public:

    // Constructor
    ExactU();

    // Destructor
    ~ExactU();

    // f(x) (must set T before using this)
    double evaluate(double x, double T);


};
#endif
