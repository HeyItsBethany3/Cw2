#include "Function1D.hpp"

/* This class specifies the initial condition for u(0,x) in our PDE problem */
class InitialU: public Function1D {

  public:
    // Constructor
    InitialU();

    // Destructor
    ~InitialU();

    // u(0,x)
    double evaluate(double x);

};
