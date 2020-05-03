#include "Elliptic.hpp"
#include "Function1.hpp"
#include "AbstractFunction.hpp"
#include <iostream>

int main(int argc, char* argv[]) {

  Function1 *f1 = new Function1();
  //Elliptic *PDE = new Elliptic(0,0,*f1,16, 1.8);
  Elliptic *PDE = new Elliptic(0,0,*f1,16, 1.8);

  (*PDE).FindSystem();
  //(*PDE).ShowSystem();
  (*PDE).FindUExact();
  //(*PDE).SolveWithTol(0.5);
  (*PDE).SolveWithIter(500);
  (*PDE).ShowApprox();
  (*PDE).ShowExact();
  (*PDE).ShowNorm();

  // Deallocate storage
  delete f1;
  delete PDE;

  return 0;
}
