
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
    void ShowSystem();

    // Solves system
    void SolveSystem();

    void ShowApprox(); // Show approximation to problem
    void ShowExact(); // Show exact solution
    void ShowNorm(); // Grid norm

    double GetNorm();

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
    double *mFvec; // Fvec

    double *uApprox; // U vector (solution)
    AbstractFunction* mFunction; // Function pointer to f


};

#endif
