
#ifndef ELLIPTIC
#define ELLIPTIC

#include "AbstractFunction.hpp"
#include <string>

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
    void ShowSystem();

    // Solves system
    void SolveWithIter(const int iter); // solves with number of iterations
    void SolveWithTol(const double tol); // Solves with specific tolerance - need to have already found UExact

    void ShowApprox(); // Show approximation to problem
    void ShowExact(); // Show exact solution
    void ShowNorm(); // Grid norm

    void UnconstrainedSol(); // Finds unconstrained solution

    double GetNorm();

    void PlotApproximation(std::string constraint);

    void FindUExact();


  protected:
    double alpha; // alpha
    double beta; // beta
    int n; // number of mesh points
    double h; // step-size
    int m; // n-1

    double *mNodes; // Interior points to solve for

    double *mDiag; // Diagonal
    double *mUpper; // Upper diagonal
    double *mLower; // Lower diagonal
    double *mFvec; // Fvec

    double *uApprox; // U vector (solution)
    AbstractFunction* mFunction; // Function pointer to f

    double w; // weight for SOR method

    double *uExact;
    double *uUnconstrained;

};

#endif
