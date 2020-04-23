#include "Functions.hpp"
#include "Option.hpp"
#include <iostream>
#include <cmath>

int main(int argc, char* argv[]) {

  Functions* f1 = new Functions(100, 0.1, 0.5);
  Option* option = new Option(100, 0.1, 0.5, 5, 300, *f1, 10000, 100);

  (*option).ConstructMatrix();
  //(*option).ShowMatrix();
  (*option).Approximate();

  //(*option).ShowApprox();
  //(*option).ShowExact();
  (*option).ShowError();
  (*option).ShowNorm();


  // Deallocate storage
  delete f1;
  delete option;
  return 0;
}
