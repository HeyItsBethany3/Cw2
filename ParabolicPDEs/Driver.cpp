#include "Parabolic.hpp"
#include "InitialU.hpp"
#include "ExactU.hpp"
#include "Function1D.hpp"
#include "Function2D.hpp"
#include <iostream>
#include <cmath>
#include <fstream>

// Use even n and l values
void plot(const double T, const int n, const int l) {
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

  (*PDE).SaveInitial();
  (*PDE2).SaveApprox();
  (*PDE).SaveApprox();

  system("cp ParabolicPlot.csv ../../../MATLAB/");

  delete uInit;
  delete uExact;
  delete PDE;
  delete PDE2;


}

void plotError(int start, int iter, double c) {
  int n = start;
  system("rm ParabolicError.csv");
  std::ofstream file;
  file.open("ParabolicError.csv");
  assert(file.is_open());
  double T = 1;
  double l = (pow(n,2)*T)/double(c);

  for(int i=1; i<=iter; i++) {
    InitialU *uInit = new InitialU(); // u_0 (initial condition)
    ExactU *uExact = new ExactU(); // exact u function
    Parabolic *PDE = new Parabolic(pow(M_PI,-2), T,0,1,*uInit, *uExact, n, l);
    (*PDE).constructMatrix();
    (*PDE).Approximate();


    // saves h, deltaT and approximation
    file << 1/double(n) << "," << (*PDE).GetMaxError() << "," << std::endl;

    n = n*2;
    l = (pow(n,2)*T)/double(c);

    delete uInit;
    delete uExact;
    delete PDE;
  }

  file.close();
  system("cp ParabolicError.csv ../../../MATLAB/");


}

void plot(const double T, const int n, const int l);

int main(int argc, char* argv[]) {

  /*
  InitialU *u0 = new InitialU(); // u_0 (initial condition)
  ExactU *uExact = new ExactU(); // exact u function
  Parabolic *PDE = new Parabolic(pow(M_PI,-2), 1,0,1,*u0, *uExact, 100, 10000);
  (*PDE).constructMatrix();
  //(*PDE).ShowMatrix();
  (*PDE).Approximate();
  //(*PDE).ShowApprox();
  //(*PDE).ShowExact();
  (*PDE).Norm();

  // Deallocate storage
  delete u0;
  delete uExact;
  delete PDE;
  */

  //plot(1, 10, 10);

  //plotError(4,8, 0.5);
  plotError(4,8, 1);
  return 0;
}
