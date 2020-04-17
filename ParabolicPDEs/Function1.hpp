#include "AbstractFunction.hpp"

// This function specifies the initial condition for u in our PDE problem
class Function1: public AbstractFunction {

  public:
    // Constructor
    Function1();

    // Destructor
    ~Function1();

    // f(x)
    double evaluate(double x);

};
