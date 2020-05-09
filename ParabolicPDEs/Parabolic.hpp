#ifndef PARABOLIC
#define PARABOLIC

#include "Function1D.hpp"
#include "Function2D.hpp"

/* Class for solving the Parabolic PDE problem in Q2 */

class Parabolic {

  public:
    // Constructor
    Parabolic(const double A, const double time, const double g_0,
      const double g_1, Function1D& InitialU, Function2D& ExactU,
      const int N, const int L);

    // Destructor
    ~Parabolic();

    // Construct matrix A such that A u_n+1 = u_n
    void constructMatrix();

    // Displays matrix
    void ShowMatrix();

    // Solves the system to find approximation for u(T,x)
    void Approximate();

    void ShowApprox(); // Shows approximation for u(T,x)
    void ShowExact(); // Shows exact value at T
    void ShowError(); // Show absolute errors
    void ShowNorm(); // Shows max error norm

    double GetMaxError(); // Retrieves max error norm

    // Saving to files
    void SaveInitial(); // Saves solution at t=0
    void SaveApprox(); // Saves approx and exact solution at T

  protected:

    // Problem set up
    double a; // conductivity coefficient
    double T; // final time value (to estimate u at)
    double g0; // u(t,0) boundary condition at x=0
    double g1; // u(t,1) boundary condition at x=1

    Function1D* mInitialU; // u_0 initial condition
    Function2D* mExactU; // exact u function
    double* uApprox; // Final approximation for u

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

};

#endif
