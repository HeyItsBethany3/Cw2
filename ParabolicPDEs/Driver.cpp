#include "Parabolic.hpp"
#include "InitialU.hpp"
#include "ExactU.hpp"
#include "Function1D.hpp"
#include "Function2D.hpp"
#include <iostream>
#include <cmath>

int main(int argc, char* argv[]) {

  InitialU *u0 = new InitialU(); // u_0 (initial condition)
  ExactU *uExact = new ExactU(); // exact u function
  Parabolic *PDE = new Parabolic(pow(M_PI,-2), 1,0,1,*u0, *uExact, 10, 2);
  (*PDE).constructMatrix();
  //(*PDE).ShowMatrix();
  (*PDE).Approximate();
  (*PDE).ShowApprox();
  (*PDE).ShowExact();


  // Deallocate storage
  delete u0;
  delete uExact;
  delete PDE;

  return 0;
}
