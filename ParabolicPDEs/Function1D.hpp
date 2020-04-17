#ifndef FUNCTION1D
#define FUNCTION1D

/* Class which specifies a function f(x) and its derivative f'(x) */

class Function1D {
  public:

    // f(x)
    virtual double evaluate(double x) = 0;

    // Destructor
    virtual ~Function1D() = 0;
};

#endif
