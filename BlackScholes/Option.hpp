#ifndef OPTIONS
#define OPTIONS

#include "AbstractFunctions.hpp"

class Option {

  public:
    // Constructor
    Option(const double strike, const double interest, const double sigma,
      const double maturity, const double maxX,
      AbstractFunctions& aFunction, const int N, const int L);

    // Destructor
    ~Option();

    // Construct matrix A
    void ConstructMatrix();

    // Show matrix
    void ShowMatrix();

    // Find approximation for u
    void Approximate();


    void ShowApprox(); // Shows approximation
    void ShowExact(); // Show exact values
    void ShowError(); // Show errors
    void Norm(); // Shows grid function norm



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

    double* uApprox; // Final approximation for u



};

#endif
