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

    // u(0,t) (Boundary function at x=0)
    double f0(double t);

    // u(R,t) (Boundary function at x=R)
    double fR(double R, double t);

};

#endif
