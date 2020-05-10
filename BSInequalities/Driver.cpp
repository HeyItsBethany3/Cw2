#include "Functions.hpp"
#include "Option.hpp"
#include <iostream>
#include <cmath>
#include <fstream>

/* Creates a plot of approximation and exact solution using maturity T,
interest rate r, stock volatilty sigma, spatial step size 1/n and time step size T/l */
void plot(const double T, const int n, const int l) {
  system("rm BSIneqPlot.csv");

  Functions* f1 = new Functions(100, 0.05, 0.5);

  // Approximates u at T
  Option* option1 = new Option(100, 0.05, 0.5, T, 300, 1.8, *f1, n, l);
  (*option1).ConstructMatrix();
  //(*option1).SolveWithIter(100);
  (*option1).SolveConvergence(0.00001);
  (*option1).FindEuropean();

  // Approximates u at T/2
  Option* option2 = new Option(100, 0.05, 0.5, T/double(2), 300, 1.8, *f1, n, l/double(2));
  (*option2).ConstructMatrix();
  //(*option2).SolveWithIter(100);
  (*option2).SolveConvergence(0.00001);
  (*option2).FindEuropean();

  // Saves data to a file
  (*option1).SaveInitial();
  (*option2).SaveApprox();
  (*option1).SaveApprox();

  system("cp BSIneqPlot.csv ../../../MATLAB/");

  delete f1;
  delete option1;
  delete option2;
}

// Plots the free boundary/ stopping time
void plotFB(const double T, const int n, const int l, const std::string param) {
  system("rm BSIneqFB.csv");

  Functions* f1 = new Functions(100, 0.05, 0.5);
  Option* option1 = new Option(100, 0.05, 0.5, T, 300, 1.8, *f1, n, l);
  (*option1).ConstructMatrix();
  //(*option1).SolveWithIter(100);
  (*option1).SolveConvergence(0.00001);
  (*option1).SaveFB(param);

  system("cp BSIneqFB.csv ../../../MATLAB/");
  delete f1;
  delete option1;
}

// Function prototypes
void plot(const double T, const int n, const int l);
void plotFB(const double T, const int n, const int l, const std::string param);

int main(int argc, char* argv[]) {

  //plot(5,100,100);
  plotFB(5,1000,300,"x");

  return 0;
}
