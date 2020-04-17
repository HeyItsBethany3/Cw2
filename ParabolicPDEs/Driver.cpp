#include "Parabolic.hpp"
#include "Function1.hpp"
#include "Function2.hpp"
#include "AbstractFunction.hpp"
#include <iostream>
#include <cmath>

int main(int argc, char* argv[]) {

  Function1 *u0 = new Function1(); // u_0 (initial condition)
  Function2 *uExact = new Function2(); // exact u function
  Parabolic *PDE = new Parabolic(pow(M_PI,-2), 1,0,1,*u0, *uExact, 5, 2);
  (*PDE).constructMatrix();
  //(*PDE).ShowMatrix();


  // Deallocate storage
  delete u0;
  delete uExact;
  delete PDE;

  return 0;
}
