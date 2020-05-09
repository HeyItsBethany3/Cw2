#ifndef ABSTRACTFUNCTIONS
#define ABSTRACTFUNCTIONS

/* Class which specifies functions for problem set up for solving the price
of an American option */

class AbstractFunctions {
  public:

    // Payoff of option given stock price at terminal time
    virtual double payoff(double x) = 0;

    // Exact value for u(x,t) where x is the stock price and t is time to maturity
    virtual double exactU(double x, double t) = 0;

    // u(0,t) (Boundary function at x=0)
    virtual double f0(double t) = 0;

    // u(R,t) (Boundary function at x=R)
    virtual double fR(double R, double t) = 0;

    // Destructor
    virtual ~AbstractFunctions() = 0;

  protected:
    double K; // Strike price
    double r; // Interest rate
    double vol; // Volatility of stock price

};

#endif
