#ifndef ABSTRACTFUNCTIONS
#define ABSTRACTFUNCTIONS

/* Class which specifies functions for problem set up for solving the price
of a European option */

class AbstractFunctions {
  public:

    // Payoff of option given stock price at terminal time
    virtual double payoff(double x) = 0;

    // Exact value for u(x,t) where x is the stock price and t is time to maturity
    virtual double exactU(double x, double t) = 0; // f_R

    // f_0 = u(t,0) (Boundary function at x=0)
    virtual double f0(double t) = 0;

    // Destructor
    virtual ~AbstractFunctions() = 0;

  protected:
    double K; // strike price
    double r; // interest rate
    double vol; // volatility of stock price

};

#endif
