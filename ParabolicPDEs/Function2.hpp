#include "AbstractFunction.hpp"

// This function specifies the exact solution for u in our PDE problem
class Function2: public AbstractFunction {

  public:
    // Constructor
    Function2();

    // Destructor
    ~Function2();

    // f(x) (must set T before using this)
    double evaluate(double x);

    // Change time step
    void changeT(const double t);

  protected:
    double T; // time to evaluate function at


};
