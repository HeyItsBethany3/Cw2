#ifndef ABSTRACTFUNCTION
#define ABSTRACTFUNCTION

/* Class which specifies a function f(x) and its derivative f'(x) */

class AbstractFunction {
  public:

    // f(x)
    virtual double evaluate(double x) = 0;

    // Destructor
    virtual ~AbstractFunction() = 0;
};

#endif
