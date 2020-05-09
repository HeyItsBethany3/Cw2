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
    //(*PDE).SolveWithTol(0.05);
    (*PDE).SolveWithIter(10e5);

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

// Creates a table of grid error values
void tableError(int start, int iter) {
  int n = start;
  system("rm EllipticIneqTable.csv");
  std::ofstream file;
  file.open("EllipticIneqTable.csv");
  assert(file.is_open());

  double* error;
  double* h;
  error = new double[iter];
  h = new double[iter];

  for(int i=1; i<=iter; i++) {
    Function1 *fun = new Function1();
    Elliptic *PDE = new Elliptic(0,0,*fun,n,1.8);
    (*PDE).FindSystem();
    (*PDE).FindUExact();
    //(*PDE).SolveWithTol(0.001);
    //(*PDE).SolveWithIter(10000);

    // Store error and mesh spacing
    error[i-1] = (*PDE).GetNorm() ;
    h[i-1] = 1/double(n);

    n = n*2;

    delete fun;
    delete PDE;
  }
  // Saves data to file
  file << "h" << ",";
  file << "Max Error" << "," << "Error[i]/Error[i-1]" << ","<< std::endl;
  for(int i=0; i< iter; i++) {
    // saves h, deltaT and approximation
    file << h[i] << ",";
    file << error[i] << ",";

    if (i>=1) {
      // Calculates difference between errors
      file << error[i-1]/double(error[i]) << ",";
    }
    file << std::endl;
  }

  delete error;
  delete h;
  file.close();
}


// Function prototypes
void plot(const int n, std::string constraint, std::string solveMethod,
  const double parameter);
void plotError(int start, int iter);
void tableError(int start, int iter);

int main(int argc, char* argv[]) {

  //plot(16, "constrained","iter", 8);
  //plot(100, "unconstrained", "tol", 0.05);
  plotError(2,6);
  //tableError(8,6);

  return 0;
}
