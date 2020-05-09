#include "Elliptic.hpp"
#include "Function1.hpp"
#include "AbstractFunction.hpp"
#include <iostream>
#include <fstream>
#include <string>

// Plots an approximation
void plot(const int n, std::string constraint, std::string solveMethod,
  const double parameter) {
  /* n+1 is the number of spatial mesh points
  solveMethod can be "iter" or "tol" - specifies whether to solve the system
  using a specific number of iterations or until the error is less than a tolerance
  constraint can be "constrained" or "unconstrained" - "constrained" plots the
  approximation against its exact solution, "unconstrained" plots it against the
  corresponding solution for Q1 (with no inequalities)
  */
  Function1 *fun = new Function1();
  Elliptic *PDE = new Elliptic(0,0,*fun,n,1.8);
  (*PDE).FindSystem();
  (*PDE).FindUExact();
  (*PDE).UnconstrainedSol();
  if (solveMethod=="iter") {
    // Solves using number of iterations
    (*PDE).SolveWithIter(parameter);
  } else {
    // Solves until approximation converges
    (*PDE).SolveWithTol(parameter);
  }
  (*PDE).PlotApproximation(constraint);

  delete fun;
  delete PDE;

}

// Creates error convergence plot
void plotError(int start, int iter) {

  int n = start; // initial n value
  system("rm EllipticIneqError.csv");
  std::ofstream file;
  file.open("EllipticIneqError.csv");
  assert(file.is_open());

  for(int i=1; i<=iter; i++) {
    Function1 *fun = new Function1();
    Elliptic *PDE = new Elliptic(0,0,*fun,n,1.8);
    (*PDE).FindSystem();
    (*PDE).FindUExact();
    (*PDE).SolveWithTol(0.05);

    // Saves mesh size and grid error norm to file
    file << 1/double(n) << "," << (*PDE).GetNorm() << "," << std::endl;

    // Updates number of mesh points
    n = n*2;

    delete fun;
    delete PDE;
  }

  file.close();
  system("cp EllipticIneqError.csv ../../../MATLAB/");

}

// Function prototypes
void plot(const int n, std::string constraint, std::string solveMethod,
  const double parameter);
void plotError(int start, int iter);

int main(int argc, char* argv[]) {

  //plot(16, "constrained","iter", 8);
  plot(16, "unconstrained", "tol", 0.05);
  //plotError(8,6);

  return 0;
}
