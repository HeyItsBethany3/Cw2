#include "Elliptic.hpp"
#include "AbstractFunction.hpp"
#include <iostream>
#include <cmath>
#include <fstream>

// We solve for the interior points u_1,...,u_m where n=m+1 and h=1/n
Elliptic::Elliptic(const double a, const double b, AbstractFunction& aFunction, const int meshPoints) {

  alpha = a;
  beta = b;
  n = meshPoints;
  h = double(1)/double(n);
  mFunction = &aFunction;

  mNodes = new double[n-1]; // x values
  mDiag = new double[n-1]; // diagonal vector for A
  mUpper = new double[n-2]; // upper diagonal vector for A
  mLower = new double[n-2]; // lower diagonal vector for A
  mFvec = new double[n-1]; // F vector
  uApprox = new double[n-1]; // Stores u approximation

  Elliptic::Nodes(); // Constructs nodes automatically

}

// Constructs nodes
void Elliptic::Nodes() {
  //std::cout << "\nNodes: ";
  for(int i=1; i<n; i++) {
    mNodes[i-1] = i*h;
    //std::cout << mNodes[i-1] << " ";
  }
}

// Constructs A matrix and F vector
void Elliptic::FindSystem() {
  // Find diagonal elements of A
  for (int i=0; i<n-1; i++) {
    mDiag[i] = -2;
  }

  // Construct upper diagonal elements of A
  for (int i=0; i<n-2; i++) {
    mUpper[i]=1;
  }

  // Construct lower diagonal elements of A
  for (int i=0; i<n-2; i++) {
    mLower[i]=1;
  }

  // Constructs F vector
  const int m = n-1;
  double factor = pow(h,2);
  mFvec[0] = -alpha-(factor*(*mFunction).evaluateF(mNodes[0]));
  mFvec[m-1] = -beta -(factor*(*mFunction).evaluateF(mNodes[m-1]));
  for(int i=1; i<m-1; i++) {
    mFvec[i] = -(factor*(*mFunction).evaluateF(mNodes[i]));
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
void Elliptic::SolveSystem() {

  // Creates dummy vectors
  double *delta, *Gvec;
  delta = new double[n-1];
  Gvec = new double[n-1];
  for(int i=0; i<=n-2; i++) {
  delta[i] = mDiag[i];
  Gvec[i] = mFvec[i];
  }

  // Elimination stage
  for(int i=1; i<=n-2; i++)
  {
    delta[i] = delta[i] - mUpper[i-1]*(mLower[i-1]/delta[i-1]);
    Gvec[i] = Gvec[i] - Gvec[i-1]*(mLower[i-1]/delta[i-1]);
  }

  // Backsolve
  uApprox[n-2] = Gvec[n-2]/delta[n-2];
  for(int i=n-3; i>=0; i--) {
    uApprox[i] = ( Gvec[i] - mUpper[i]*uApprox[i+1] )/delta[i];
  }

  // Deallocates storage
  delete delta;
  delete Gvec;
}

// Show approximation
void Elliptic::ShowApprox() {
  std::cout << "\nApproximation: ";
  for (int i=0; i<n-1; i++) {
    std::cout << uApprox[i] << " ";
  }
  std::cout << std::endl;
}

// Show exact solution
void Elliptic::ShowExact() {
  std::cout << "\nExact solution: ";
  for (int i=0; i<n-1; i++) {
    std::cout << (*mFunction).exactU(mNodes[i]) << " ";
  }
  std::cout << std::endl;
}

// Shows the grid norm
void Elliptic::ShowNorm() {
  double sum = 0;
  for(int i=0; i<n-1; i++) {
    sum = sum +  fabs(uApprox[i]-(*mFunction).exactU(mNodes[i]));
  }
  sum = sqrt(sum *h);
  std::cout << "\nGrid norm: " << sum << "\n";
}

// Retrieves the grid norm
double Elliptic::GetNorm() {
  double sum = 0;
  for(int i=0; i<n-1; i++) {
    sum = sum +  fabs(uApprox[i]-(*mFunction).exactU(mNodes[i]));
  }
  sum = sqrt(sum *h);
  return sum;
}

void Elliptic::PlotApproximation() {
  // Writes data to a file
  std::ofstream file;
  file.open("EllipticPlot.csv");
  assert(file.is_open());

  // Saves x values
  file << 0 << ",";
  for(int i=0; i<n-1; i++) {
    file << mNodes[i] << ",";
  }
  file << 1 << ",";
  file << std::endl;

  // Saves Approximation
  file << alpha << ",";
  for(int i=0; i<n-1; i++) {
    file << uApprox[i] << ",";
  }
  file << beta << ",";
  file << std::endl;

  // Saves exact solution
  file << (*mFunction).exactU(0) << ",";
  for(int i=0; i<n-1; i++) {
    file << (*mFunction).exactU(mNodes[i]) << ",";
  }
  file << (*mFunction).exactU(1) << ",";
  file.close();
  system("cp EllipticPlot.csv ../../../MATLAB/");
}

// Deallocates storage
Elliptic::~Elliptic() {
  delete mNodes;
  delete mDiag;
  delete mUpper;
  delete mLower;
  delete mFvec;
  delete uApprox;
}
