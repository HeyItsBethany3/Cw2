#include "Parabolic.hpp"
#include "InitialU.hpp"
#include "ExactU.hpp"
#include <iostream>
#include <cmath>
#include <fstream>

/* Creates a plot of approximation and exact solution using maturity T,
spatial step size 1/n and time step size T/l */
void plot(const double T, const int n, const int l) {
  // Always use even n and l values
  system("rm ParabolicPlot.csv");

  InitialU *uInit = new InitialU(); // u_0 (initial condition)
  ExactU *uExact = new ExactU(); // exact u function

  // Approximates u at T
  Parabolic *PDE = new Parabolic(pow(M_PI,-2), T,0,1,*uInit, *uExact, n, l);
  (*PDE).constructMatrix();
  (*PDE).Approximate();

  // Approximates u at T/2
  Parabolic *PDE2 = new Parabolic(pow(M_PI,-2), T/double(2),0,1,*uInit, *uExact,
  n, l/double(2));
  (*PDE2).constructMatrix();
  (*PDE2).Approximate();

  // Saves data to a file
  (*PDE).SaveInitial();
  (*PDE2).SaveApprox();
  (*PDE).SaveApprox();

  system("cp ParabolicPlot.csv ../../../MATLAB/");

  delete uInit;
  delete uExact;
  delete PDE;
  delete PDE2;
}

// Creates an error convergence plot
void plotError(int start, int iter, double c) {
  int n = start;
  system("rm ParabolicError.csv");
  std::ofstream file;
  file.open("ParabolicError.csv");
  assert(file.is_open());
  double T = 1; // Maturity T
  double l = (pow(n,2)*T)/double(c); // Finds corresponding time step size

  for(int i=1; i<=iter; i++) {
    InitialU *uInit = new InitialU(); // u_0 (initial condition)
    ExactU *uExact = new ExactU(); // exact u function
    Parabolic *PDE = new Parabolic(pow(M_PI,-2), T,0,1,*uInit, *uExact, n, l);
    (*PDE).constructMatrix();
    (*PDE).Approximate();

    // For each n, save the h (mesh size) and max norm to file
    file << 1/double(n) << "," << (*PDE).GetMaxError() << "," << std::endl;

    // Change mesh spacing
    n = n*2;
    l = (pow(n,2)*T)/double(c);

    delete uInit;
    delete uExact;
    delete PDE;
  }
  file.close();
  system("cp ParabolicError.csv ../../../MATLAB/");
}

// Creates a table of max error values
void tableError(int start, int iter, double c) {
  int n = start;
  system("rm ParabolicTable.csv");
  std::ofstream file;
  file.open("ParabolicTable.csv");
  assert(file.is_open());
  double T = 1;  // Maturity T
  double l = (pow(n,2)*T)/double(c); // Finds corresponding time step size

  double* error;
  double* h;
  double* t;
  error = new double[iter];
  h = new double[iter];
  t = new double[iter];

  for(int i=1; i<=iter; i++) {
    InitialU *uInit = new InitialU(); // u_0 (initial condition)
    ExactU *uExact = new ExactU(); // exact u function
    Parabolic *PDE = new Parabolic(pow(M_PI,-2), T,0,1,*uInit, *uExact, n, l);
    (*PDE).constructMatrix();
    (*PDE).Approximate();

    // Store error and mesh spacing
    error[i-1] = (*PDE).GetMaxError();
    h[i-1] = 1/double(n);
    t[i-1] = T/double(l);

    // Change mesh spacing
    n = n*2;
    l = (pow(n,2)*T)/double(c);

    delete uInit;
    delete uExact;
    delete PDE;
  }

  // Saves data to file
  file << "h" << "," << "Delta t" << ",";
  file << "Max Error" << "," << "Error[i]/Error[i-1]" << ","<< std::endl;
  for(int i=0; i< iter; i++) {
    // saves h, deltaT and approximation
    file << h[i] << "," << t[i] << ",";
    file << error[i] << ",";

    if (i>=1) {
      // Calculates difference between errors
      file << error[i]/double(error[i-1]) << ",";
    }
    file << std::endl;
  }
  file.close();

  delete error;
  delete h;
  delete t;

  system("cp ParabolicTable.csv ../../../MATLAB/");
}


// Function Prototypes
void plot(const double T, const int n, const int l);
void plotError(int start, int iter, double c);
void tableError(int start, int iter, double c);

int main(int argc, char* argv[]) {

  //plot(1, 10, 10);
  //plotError(4,8, 1);
  tableError(4,8,1);

  return 0;
}
