#include "AbstractFunction.hpp"

class Function1: public AbstractFunction {

  public:
    // Constructor
    Function1();

    // Destructor
    ~Function1();

    // f(x)
    double f(double x);

    // psi(x)
    double psi(double x);

    // Initial guess for u at all interior points
    double init(double x);

    // exact value for u(x)
    double exactU(double x);

};
