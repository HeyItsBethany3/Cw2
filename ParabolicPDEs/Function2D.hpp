#ifndef FUNCTION2D
#define FUNCTION2D

/* Abstract class for function with two inputs */

class Function2D {
  public:

    // f(x)
    virtual double evaluate(double x, double T) = 0;

    // Destructor
    virtual ~Function2D() = 0;
};

#endif
