
#ifndef ELLIPTIC
#define ELLIPTIC

#include "AbstractFunction.hpp"
#include <string>

/* Class to solve the elliptic inequality PDE problem */

class Elliptic {

  public:

    // Constructor
    Elliptic(const double a, const double b, AbstractFunction& aFunction,
      const int meshPoints, const double weight);

    // Destructor
    ~Elliptic();

    // Constructs nodes
    void Nodes();

    // Constructs system
    void FindSystem();

    // Show system to solve
    void ShowSystem();

    // Solves system
    void SolveWithIter(const int iter); // Solves with number of iterations
    void SolveWithTol(const double tol); // Solves with specific tolerance

    void ShowApprox(); // Show approximation to problem
    void ShowExact(); // Show exact solution
    void ShowNorm(); // Show grid norm

    double GetNorm(); // Retrieves grid norm

    void UnconstrainedSol(); // Finds unconstrained solution (for Q1)
    void FindUExact(); // Finds exact solution

    // Saves data to file to create plot
    void PlotApproximation(std::string constraint);


  protected:
    double alpha; // alpha
    double beta; // beta
    int n; // n+1 is the number of mesh points
    double h; // spatial step-size
    int m; // n-1

    double *mNodes; // Interior points to solve for

    double *mDiag; // Diagonal
    double *mUpper; // Upper diagonal
    double *mLower; // Lower diagonal
    double *mFvec; // Fvec

    double *uApprox; // Approximation for u
    AbstractFunction* mFunction; // Function pointer to f

    double w; // weight for SOR method

    double *uExact; // Exact u(x)
    double *uUnconstrained; // Unconstrained solution

};

#endif
