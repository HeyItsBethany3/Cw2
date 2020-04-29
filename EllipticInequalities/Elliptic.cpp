#include "Elliptic.hpp"
#include "AbstractFunction.hpp"
#include <iostream>
#include <cmath>

// n=m+1 mesh points, h=1/n
Elliptic::Elliptic(const double a, const double b, AbstractFunction& aFunction,
  const int meshPoints, const double weight) {
  alpha = a;
  beta = b;
  n = meshPoints;
  m = n-1;
  h = double(1)/double(n);
  w = weight;
  mFunction = &aFunction;

  mNodes = new double[n-1];
  mDiag = new double[n-1]; // diagonal vector
  mUpper = new double[n-2]; // upper diagonal vector
  mLower = new double[n-2]; // lower diagonal vector
  mFvec = new double[n-1]; // F vector
  uApprox = new double[n-1]; // U solution
  uExact = new double[n-1];

  Elliptic::Nodes(); // Constructs nodes automatically

}

void Elliptic::Nodes() {
  for(int i=1; i<n; i++) {
    mNodes[i-1] = i*h;
  }
}

void Elliptic::FindSystem() {
  // Find diagonal elements of A
  for (int i=0; i<n-1; i++) {
    mDiag[i] = 2;
  }

  // Construct upper diagonal elements of A
  for (int i=0; i<n-2; i++) {
    mUpper[i]=-1;
  }

  // Construct lower diagonal elements of A
  for (int i=0; i<n-2; i++) {
    mLower[i]=-1;
  }

  // Constructs F vector
  const int m = n-1;
  double factor = pow(h,2);
  mFvec[0] = alpha+(factor*(*mFunction).f(mNodes[0]));
  mFvec[m-1] = beta +(factor*(*mFunction).f(mNodes[m-1]));
  for(int i=1; i<m-1; i++) {
    mFvec[i] = (factor*(*mFunction).f(mNodes[i]));
  }

}

// Displays system to solve
void Elliptic::ShowSystem() {
  std::cout << "\nF(x): ";
  for (int i=0; i<=n-2; i++) {
    std::cout << mFvec[i] << " ";
  }

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

// Solves system of equations
void Elliptic::SolveSystem(const int iter) {

  double *uApproxOld;
  uApproxOld = new double[m];
  for(int i=0; i<m; i++) {
    uApproxOld[i] = (*mFunction).init(mNodes[i]);
  }

  double uVal, psi, psiVal;

  for(int k=1; k<=iter; k++) {
    for(int i=0; i<m; i++) {

      if (i==0) {
        uVal = (mFvec[0]-(mUpper[0]*uApproxOld[1]))/mDiag[0];
      } else if (i==m-1) {
        uVal = (mFvec[m-1]-(mLower[m-2]*uApprox[m-2]))/mDiag[m-1];
      } else {
        uVal = (mFvec[i]-(mLower[i-1]*uApprox[i-1])-(mUpper[i]*uApproxOld[i+1]))/mDiag[i];
      }
      psi = (*mFunction).psi(mNodes[i]);
      psiVal = (w*uVal)+((1-w)*uApproxOld[i]);
      if (psi<psiVal) {
        uApprox[i] = psi;
      } else {
        uApprox[i] = psiVal;
      }
    }

    for (int i=0; i<m; i++) {
      uApproxOld[i] = uApprox[i];
    }

  }

  delete uApproxOld;
}

// Solves AU=F
void Elliptic::FindUExact() {

  double* uArray;
  uArray = new double[n-1];

  //Create delta and Gvec vectors of the Triangular system
  double *delta, *Gvec;
  delta = new double[n-1];
  Gvec = new double[n-1];
  for(int i=0; i<=n-2; i++)
  {
  delta[i] = mDiag[i];
  Gvec[i] = mFvec[i];
  }

  // Elimination stage
  for(int i=1; i<=n-2; i++)
  {
    delta[i] = delta[i] - mUpper[i-1]*(mLower[i-1]/delta[i-1]);
    Gvec[i] = Gvec[i] - Gvec[i-1]*(mLower[i-1]/delta[i-1]);
  }

  //Backsolve
  uArray[n-2] = Gvec[n-2]/delta[n-2];
  for(int i=n-3; i>=0; i--)
  {
    uArray[i] = ( Gvec[i] - mUpper[i]*uApprox[i+1] )/delta[i];
  }

  for(int i=0; i<n-1; i++){
    // Free boundaries
    double x1 = sqrt(3)/double(5);
    double x2 = 1-(sqrt(3)/double(5));
    if (mNodes[i] < x1) {
      uExact[i] = uArray[i];
    } else if (mNodes[i] > x2) {
      uExact[i] = uArray[i];
    } else {
      uExact[i] = (*mFunction).psi(mNodes[i]);
    }

  }

  // Deallocates storage
  delete delta;
  delete Gvec;
  delete uArray;
}

void Elliptic::ShowApprox() {
  std::cout << "\nApproximation: ";
  for (int i=0; i<n-1; i++) {
    std::cout << uApprox[i] << " ";
  }
  std::cout << std::endl;
}

void Elliptic::ShowExact() {
  std::cout << "\nExact solution: ";
  for (int i=0; i<n-1; i++) {
    std::cout << uExact[i] << " ";
  }
  std::cout << std::endl;
}

// Shows the grid norm
void Elliptic::Norm() {
  double sum = 0;
  for(int i=0; i<n-1; i++) {
    sum = sum +  fabs(uApprox[i]-uExact[i]);
  }
  sum = sqrt(sum *h);
  std::cout << "\nGrid norm: " << sum << "\n";

}

Elliptic::~Elliptic() {
  delete mNodes;
  delete mDiag;
  delete mUpper;
  delete mLower;
  delete mFvec;
  delete uApprox;
}
