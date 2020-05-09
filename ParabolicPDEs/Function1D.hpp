#ifndef FUNCTION1D
#define FUNCTION1D

/* Abstract class for function with one input */

class Function1D {
  public:

    // f(x)
    virtual double evaluate(double x) = 0;

    // Destructor
    virtual ~Function1D() = 0;
};

#endif
