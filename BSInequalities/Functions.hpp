#ifndef FUNCTIONS
#define FUNCTIONS
#include "AbstractFunctions.hpp"

class Functions: public AbstractFunctions {

  public:
    // Constructor
    Functions(const double strike, const double interest, const double sigma);

    // Destructor
    ~Functions();

    // Payoff of option given stock price at terminal time
    double payoff(double x);

    // Exact value for u(x,t) where x is the stock price and t is time to maturity
    double exactU(double x, double t); // f_R

    // f_0 = u(t,0) (Boundary function at x=0)
    double f0(double t);
    double fR(double t);

};

#endif
