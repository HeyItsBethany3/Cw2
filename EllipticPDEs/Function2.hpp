#include "AbstractFunction.hpp"

class Function2: public AbstractFunction {

  public:
    // Constructor
    Function2();

    // Destructor
    ~Function2();

    // f(x)
    double evaluateF(double x);

    // Exact value for u(x)
    double exactU(double x);

};
