#include "Elliptic.hpp"
#include "Function1.hpp"
#include "AbstractFunction.hpp"
#include <iostream>
#include <fstream>
#include <string>


void plot(const int n, std::string constraint) {
  Function1 *fun = new Function1();
  Elliptic *PDE = new Elliptic(0,0,*fun,n,1.8);
  (*PDE).FindSystem();
  (*PDE).FindUExact();
  (*PDE).UnconstrainedSol();
  (*PDE).SolveWithIter(500); // CHANGE
  (*PDE).PlotApproximation(constraint);
  //(*PDE).ShowApprox();
  //(*PDE).ShowExact();

  delete fun;
  delete PDE;

}

void plotError(int start, int iter) {

  int n = start;
  system("rm EllipticIneqError.csv");
  std::ofstream file;
  file.open("EllipticIneqError.csv");
  assert(file.is_open());

  for(int i=1; i<=iter; i++) {
    Function1 *fun = new Function1();
    Elliptic *PDE = new Elliptic(0,0,*fun,n,1.8);
    (*PDE).FindSystem();
    (*PDE).FindUExact();
    (*PDE).SolveWithIter(100); // CHANGE



    file << 1/double(n) << "," << (*PDE).GetNorm() << "," << std::endl;

    n = n*2;

    delete fun;
    delete PDE;
  }

  file.close();
  system("cp EllipticIneqError.csv ../../../MATLAB/");

}


int main(int argc, char* argv[]) {
  /*
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

  delete f1;
  delete PDE;
  */
  // Plot an approximation
  plot(100, "constrained");

  // Plot errors
  //plotError(2,20);


  return 0;
}
