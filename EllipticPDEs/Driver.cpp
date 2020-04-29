#include "Elliptic.hpp"
#include "Function1.hpp"
#include "AbstractFunction.hpp"
#include <iostream>
#include <fstream>

void plot(const int n) {
  Function1 *fun = new Function1();
  Elliptic *pde1 = new Elliptic(0,1,*fun,n);
  (*pde1).FindSystem();
  (*pde1).SolveSystem();
  (*pde1).PlotApproximation();

  delete fun;
  delete pde1;

}

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

    file << 1/double(n) << "," << (*pde1).GetNorm() << "," << std::endl;

    n = n*2;

    delete fun;
    delete pde1;
  }

  file.close();
  system("cp EllipticError.csv ../../../MATLAB/");


}


void plot(const int n);
void plotError(int start, int iter);

int main(int argc, char* argv[]) {

  /*
  Function1 *f1 = new Function1();
  Elliptic *PDE = new Elliptic(0,1,*f1,4);
  (*PDE).FindSystem();
  //(*PDE).ShowSystem();
  (*PDE).SolveSystem();
  (*PDE).ShowApprox();
  (*PDE).ShowExact();
  (*PDE).ShowNorm();
  delete f1;
  delete PDE;

  */

  // Plot an approximation
  plot(7);

  // Plot errors
  plotError(2,20);






  return 0;
}
