#include "Functions.hpp"
#include "Option.hpp"
#include <iostream>
#include <cmath>
#include <fstream>


// Use even n and l values
void plot(const double T, const int n, const int l) {
  system("rm BSPlot.csv");

  Functions* f1 = new Functions(100, 0.1, 0.5);

  // Approximates u at T
  Option* option1 = new Option(100, 0.1, 0.5, 5, T, *f1, n, l);
  (*option1).ConstructMatrix();
  (*option1).Approximate();


  // Approximates u at T/2
  Option* option2 = new Option(100, 0.1, 0.5, 5, T/double(2), *f1, n, l/double(2));
  (*option2).ConstructMatrix();
  (*option2).Approximate();


  // Saves to file
  (*option1).SaveInitial();
  (*option2).SaveApprox();
  (*option1).SaveApprox();

  system("cp BSPlot.csv ../../../MATLAB/");

  delete f1;
  delete option1;
  delete option2;

}

void plotError(int start, int iter, double c) {
  int n = start;
  system("rm BSError.csv");
  std::ofstream file;
  file.open("BSError.csv");
  assert(file.is_open());
  double T = 1;
  double l = (pow(n,2)*T)/double(c);

  for(int i=1; i<=iter; i++) {
    Functions* f1 = new Functions(100, 0.1, 0.5);
    Option* option1 = new Option(100, 0.1, 0.5, 5, T, *f1, n, l);
    (*option1).ConstructMatrix();
    (*option1).Approximate();


    // saves h, deltaT and approximation
    file << 1/double(n) << "," << (*option1).GetMaxError() << "," << std::endl;

    n = n*2;
    l = (pow(n,2)*T)/double(c);

    delete f1;
    delete option1;
  }

  file.close();
  system("cp BSError.csv ../../../MATLAB/");


}

void plot(const double T, const int n, const int l);


int main(int argc, char* argv[]) {
  /*
  Functions* f1 = new Functions(100, 0.1, 0.5);
  Option* option = new Option(100, 0.1, 0.5, 5, 300, *f1, 10000, 100);

  (*option).ConstructMatrix();
  //(*option).ShowMatrix();
  (*option).Approximate();

  //(*option).ShowApprox();
  //(*option).ShowExact();
  (*option).ShowError();
  (*option).ShowNorm();
  delete f1;
  delete option;
  */

  //plot(10,10,300);
  plotError(4,8,1);



  return 0;
}
