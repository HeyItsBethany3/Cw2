#include "Elliptic.hpp"
#include "Function1.hpp"
#include "Function2.hpp"
#include "AbstractFunction.hpp"
#include <iostream>
#include <fstream>

// Creates a plot of approximation and exact solution
void plot(const int n) {
  Function1 *fun = new Function1();
  Elliptic *pde1 = new Elliptic(0,1,*fun,n);
  (*pde1).FindSystem();
  (*pde1).SolveSystem();
  (*pde1).PlotApproximation();

  delete fun;
  delete pde1;

}

// Creates an error convergence plot
void plotError(int start, int iter) {
  int n = start;
  system("rm EllipticError.csv");
  std::ofstream file;
  file.open("EllipticError.csv");
  assert(file.is_open());

  for(int i=1; i<=iter; i++) {
    Function1 *fun = new Function1();
    Elliptic *pde1 = new Elliptic(0,1,*fun,n);
    (*pde1).FindSystem();
    (*pde1).SolveSystem();

    // For each n, save the h (mesh size) and grid norm to file
    file << 1/double(n) << "," << (*pde1).GetNorm() << "," << std::endl;

    n = n*2;

    delete fun;
    delete pde1;
  }

  file.close();
  system("cp EllipticError.csv ../../../MATLAB/");
}

// Creates a table of grid error values
void tableError(int start, int iter) {
  int n = start;
  system("rm EllipticTable.csv");
  std::ofstream file;
  file.open("EllipticTable.csv");
  assert(file.is_open());

  double* error;
  double* h;
  error = new double[iter];
  h = new double[iter];

  for(int i=1; i<=iter; i++) {
    Function1 *fun = new Function1();
    Elliptic *pde1 = new Elliptic(0,1,*fun,n);
    (*pde1).FindSystem();
    (*pde1).SolveSystem();

    // Store error and mesh spacing
    error[i-1] = (*pde1).GetNorm() ;
    h[i-1] = 1/double(n);

    n = n*2;

    delete fun;
    delete pde1;
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
void plot(const int n);
void plotError(int start, int iter);
void tableError(int start, int iter);

int main(int argc, char* argv[]) {

  //plot(5);
  //plotError(2,15);
  tableError(2,15);

  return 0;
}
