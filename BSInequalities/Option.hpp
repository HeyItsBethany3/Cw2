#ifndef OPTIONS
#define OPTIONS

#include "AbstractFunctions.hpp"
#include <string>

/* Class for finding the price of a American option in Q5 */

class Option {

  public:
    // Constructor
    Option(const double strike, const double interest, const double sigma,
          const double maturity, const double maxX, const double weight,
          AbstractFunctions& aFunction, const int N, const int L);

    // Destructor
    ~Option();

    // Construct matrix A
    void ConstructMatrix();

    // Show matrix
    void ShowMatrix();

    // Find corresponding European approximation (Q3)
    void FindEuropean();

    // Find approximation for u(x,T)
    void SolveWithIter(const int iter);

    // Solve until approximation converges
    void SolveConvergence(const double tol);

    void ShowApprox(); // Shows approximation for u(T,x)
    void ShowExact(); // Shows exact value at T
    void ShowError(); // Show absolute errors
    void ShowNorm(); // Shows max error norm

    double GetMaxError(); // Retrieves max error norm

    // Save to files
    void SaveInitial(); // Saves solution at t=0
    void SaveApprox(); // Saves approx and exact solution at T
    void SaveFB(const std::string parameter); // Saves free boundary


  protected:
    double K; // strike price
    double r; // interest rate
    double vol; // volatility
    double T; // length of maturity
    double R; // max stock price
    double w;

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

    double* uApprox; // Final approximation for u (American option)
    double* European; // Stores European option value
    double* FreeBoundaryT; // Stores stopping times for each t
    double* FreeBoundaryX;
    double* FBNotFoundT;
    double* FBNotFoundX;

    //double* FBNotFound; // Checks whether a stopping time has been found yet


};

#endif
