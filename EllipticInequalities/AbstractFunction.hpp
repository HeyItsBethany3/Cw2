#ifndef ABSTRACTFUNCTION
#define ABSTRACTFUNCTION

/* Class which specifies a function f(x) and its derivative f'(x) */

class AbstractFunction {
  public:

    // f(x)
    virtual double f(double x) = 0;

    // Psi(x)
    virtual double psi(double x) = 0;

    // Initial guess for u at all interior points
    virtual double init(double x) = 0;

    // Destructor
    virtual ~AbstractFunction() = 0;
};

#endif