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

  // Approximate U at 0 (Initial condition)



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

  plot(1, 10, 10);
  return 0;
}
