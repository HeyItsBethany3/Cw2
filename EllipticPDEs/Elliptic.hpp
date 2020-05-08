
#ifndef ELLIPTIC
#define ELLIPTIC

#include "AbstractFunction.hpp"

class Elliptic {

  public:

    // Constructor
    Elliptic(const double a, const double b, AbstractFunction& aFunction, const int meshPoints);

    // Destructor
    ~Elliptic();

    // Constructs nodes
    void Nodes();

    // Constructs system
    void FindSystem();
    // Displays system to solve
    void ShowSystem();

    // Solves system (finds approximation)
    void SolveSystem();

    void ShowApprox(); // Show approximation
    void ShowExact(); // Show exact solution
    void ShowNorm(); // Shows grid norm

    double GetNorm(); // Retrieves grid norm

    // Writes approximation and exact u to file
    void PlotApproximation();

  protected:
    double alpha; // alpha
    double beta; // beta
    int n; // number of mesh points
    double h; // step-size

    double *mNodes; // Interior points to solve for

    double *mDiag; // Diagonal
    double *mUpper; // Upper diagonal
    double *mLower; // Lower diagonal
    double *mFvec;

    double *uApprox; // Approximation for u
    AbstractFunction* mFunction; // Function pointer to f


};

#endif
