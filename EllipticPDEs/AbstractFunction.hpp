#ifndef ABSTRACTFUNCTION
#define ABSTRACTFUNCTION

/* Class which specifies f(x) and the exact solution u(x) */

class AbstractFunction {
  public:

    // f(x)
    virtual double evaluateF(double x) = 0;

    // Exact value for u(x)
    virtual double exactU(double x) = 0;

    // Destructor
    virtual ~AbstractFunction() = 0;
};

#endif
