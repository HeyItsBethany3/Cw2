#include "Functions.hpp"
#include "Option.hpp"
#include <iostream>
#include <cmath>
#include <fstream>

/* Creates a plot of approximation and exact solution using maturity T,
interest rate r, stock volatilty sigma, spatial step size 1/n and time step size T/l */
void plot(const double T, const double r, const double sigma, const int n, const int l) {
  system("rm BSPlot.csv");

  Functions* f1 = new Functions(100, r, sigma);

  // Approximates u at T
  Option* option1 = new Option(100, r, sigma, T, 300, *f1, n, l);
  (*option1).ConstructMatrix();
  (*option1).Approximate();

  // Approximates u at T/2
  Option* option2 = new Option(100, r, sigma, T/double(2), 300, *f1, n, l/double(2));
  (*option2).ConstructMatrix();
  (*option2).Approximate();

  // Saves data to a file
  (*option1).SaveInitial();
  (*option2).SaveApprox();
  (*option1).SaveApprox();

  system("cp BSPlot.csv ../../../MATLAB/");

  delete f1;
  delete option1;
  delete option2;
}

// Creates an error convergence plot
void plotError(const double r, const double sigma, const int start,
  const int iter, const double c) {
  // 'Start' specifies n start value and 'iter' the number of results to obtain
  int n = start;
  system("rm BSError.csv");
  std::ofstream file;
  file.open("BSError.csv");
  assert(file.is_open());
  double T = 5; // Maturity T
  double l = (pow(n,1)*T)/double(c); // Finds corresponding time step size

  for(int i=1; i<=iter; i++) {
    Functions* f1 = new Functions(100, r, sigma);
    Option* option1 = new Option(100, r, sigma, T, 300, *f1, n, l);
    (*option1).ConstructMatrix();
    (*option1).Approximate();

    // For each n, save the h (mesh size) and max norm to file
    file << 1/double(n) << "," << (*option1).GetMaxError() << "," << std::endl;

    // Change mesh spacing
    n = n*2;
    l = (pow(n,1)*T)/double(c);

    delete f1;
    delete option1;
  }
  file.close();
  system("cp BSError.csv ../../../MATLAB/");
}

// Creates a table of max error values
void tableError(const double r, const double sigma, const int start,
  const int iter, const double c) {
  // 'Start' specifies n start value and 'iter' the number of results to obtain
  int n = start;
  system("rm BSTable.csv");
  std::ofstream file;
  file.open("BSTable.csv");
  assert(file.is_open());
  double T = 5; // Maturity T
  double l = (pow(n,2)*T)/double(c); // Finds corresponding time step size

  double* error;
  double* h;
  double* t;
  error = new double[iter];
  h = new double[iter];
  t = new double[iter];

  for(int i=1; i<=iter; i++) {
    Functions* f1 = new Functions(100, r, sigma);
    Option* option1 = new Option(100, r, sigma, T, 300, *f1, n, l);
    (*option1).ConstructMatrix();
    (*option1).Approximate();

    // Store error and mesh spacing
    error[i-1] = (*option1).GetMaxError();
    h[i-1] = 1/double(n);
    t[i-1] = T/double(l);

    // Change mesh spacing
    n = n*2;
    l = (pow(n,2)*T)/double(c);

    delete f1;
    delete option1;
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
      file << error[i-1]/double(error[i]) << ",";
    }
    file << std::endl;
  }
  file.close();

  delete error;
  delete h;
  delete t;
  file.close();
}

// Function prototypes
void plot(const double T, const double r, const double sigma, const int n, const int l);
void plotError(const double r, const double sigma, const int start,
  const int iter, const double c);
void tableError(const double r, const double sigma, const int start,
  const int iter, const double c);

int main(int argc, char* argv[]) {

  //plot(5,0,0.5,10,10);
  //tableError(0,0.5,8,7,1);
  //plotError(0.1,0.1,64,6,1);
  tableError(0.1,0.1,8,7,1);

  return 0;
}
