#include "Parabolic.hpp"
#include <iostream>
#include <cmath>

// Constructor
Parabolic::Parabolic(const double A, const double time, const double g_0,
  const double g_1, AbstractFunction& InitialU, AbstractFunction& ExactU,
  const int N, const int L) {
    a = A;
    T = time;
    g0 = g_0;
    g1 = g_1;
    mInitialU = &InitialU;
    mExactU = &ExactU;
    n=N;
    m=N-1;
    l=L;
    deltaT = T/double(l);

    xNodes = new double[m];
    h = 1/double(N); // spatial step-size
    for(int i=0; i<m; i++) {
      xNodes[i]=(i+1)*h;
    }

    mDiag = new double[m];
    mUpper = new double[m-1];
    mLower = new double[m-1];



  }

void Parabolic::constructMatrix() {
const double lamda = double(a*deltaT)/double(pow(h,2));

  // Diagonal elements of A
  for (int i=0; i<m; i++) {
    mDiag[i] = (lamda*2)+1;
  }

  // Construct upper diagonal elements of A
  for (int i=0; i<m-1; i++) {
    mUpper[i]=-1*lamda;
  }

  // Construct lower diagonal elements of A
  for (int i=0; i<m-1; i++) {
    mLower[i]=-1*lamda;
  }
}

// Displays system to solve
void Parabolic::ShowMatrix() {

  std::cout << "\nd: ";
  for (int i=0; i<=n-2; i++) {
    std::cout << mDiag[i] << " ";
  }

  std::cout << "\nu: ";
  for (int i=0; i<n-2; i++) {
    std::cout << mUpper[i] << " ";
  }

  std::cout << "\nl: ";
  for (int i=0; i<n-2; i++) {
    std::cout << mLower[i] << " ";
  }
  std::cout << std::endl;
}

// Destructor
Parabolic::~Parabolic() {
  delete xNodes;
  delete mDiag;
  delete mUpper;
  delete mLower;
  delete uApprox;
}
