#include "Function1D.hpp"

// This function specifies the initial condition for u in our PDE problem
class InitialU: public Function1D {

  public:
    // Constructor
    InitialU();

    // Destructor
    ~InitialU();

    // f(x)
    double evaluate(double x);

};
