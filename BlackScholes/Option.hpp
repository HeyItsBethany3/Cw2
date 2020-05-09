#ifndef OPTIONS
#define OPTIONS

#include "AbstractFunctions.hpp"

/* Class for finding the price of a European option in Q3 */

class Option {

  public:
    // Constructor
    Option(const double strike, const double interest, const double sigma,
      const double maturity, const double maxX,
      AbstractFunctions& aFunction, const int N, const int L);

    // Destructor
    ~Option();

    // Construct matrix A such that A u_n+1 = u_n + f
    void ConstructMatrix();
    // Show matrix A
    void ShowMatrix();

    // Find approximation for u(x,T)
    void Approximate();

    void ShowApprox(); // Shows approximation for u(T,x)
    void ShowExact(); // Shows exact value at T
    void ShowError(); // Show absolute errors
    void ShowNorm(); // Shows max error norm

    double GetMaxError(); // Retrieves max error norm

    // Save to files
    void SaveInitial(); // Saves solution at t=0
    void SaveApprox(); // Saves approx and exact solution at T


  protected:
    double K; // strike price
    double r; // interest rate
    double vol; // volatility
    double T; // length of maturity
    double R; // max stock price

    AbstractFunctions* mFunction; // object with useful functions

    // Discretisation choices
    int n; // number of spatial mesh points
    int m; // m=n-1
    int l; // number of time mesh points
    double h; // spatial step size
    double deltaT; // time step size

    // Nodes
    double* xNodes; // interior spatial nodes x1,...,xm

    // For matrix A
    double *mDiag; // Diagonal
    double *mUpper; // Upper diagonal
    double *mLower; // Lower diagonal

    double* uApprox; // Final approximation for u(x,T)
};

#endif
