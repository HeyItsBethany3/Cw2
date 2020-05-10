#include "Elliptic.hpp"
#include "AbstractFunction.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

// We solve for the interior points u_1,...,u_m where n=m+1 and h=1/n
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
  mDiag = new double[n-1]; // diagonal vector for A
  mUpper = new double[n-2]; // upper diagonal vector for A
  mLower = new double[n-2]; // lower diagonal vector for A
  mFvec = new double[n-1]; // F vector
  uApprox = new double[n-1]; // Stores u approximation
  uExact = new double[n-1]; // Stores exact u(x)
  uUnconstrained = new double[n-1]; // Stores unconstrained approximation

  Elliptic::Nodes(); // Constructs nodes automatically

}

void Elliptic::Nodes() {
  for(int i=1; i<n; i++) {
    mNodes[i-1] = i*h;
  }
}

// Constructs A matrix and F vector
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

// Solves problem using "iter" iterations of the SOR method
void Elliptic::SolveWithIter(const int iter) {

  double *uApproxOld; // Stores previous iteration of method
  uApproxOld = new double[m];
  for(int i=0; i<m; i++) {
    uApproxOld[i] = (*mFunction).init(mNodes[i]);
    uApprox[i] = (*mFunction).init(mNodes[i]);
  }

  double uVal, psi, psiVal;

  for(int k=1; k<=iter; k++) {

    // Implements SOR method
    for(int i=0; i<m; i++) {

      if (i==0) {
        uVal = (mFvec[0]-(mUpper[0]*uApproxOld[1]))/mDiag[0];
      } else if (i==m-1) {
        uVal = (mFvec[m-1]-(mLower[m-2]*uApprox[m-2]))/mDiag[m-1];
      } else {
        uVal = (mFvec[i]-(mLower[i-1]*uApprox[i-1])-(mUpper[i]*uApproxOld[i+1]))/mDiag[i];
      }
      psi = (*mFunction).psi(mNodes[i]);

      psiVal = (w*uVal)+((1.0-w)*uApproxOld[i]);

      // Updates u value
      if ( psi <psiVal) {
        uApprox[i] = psi;
      } else {
        uApprox[i] = psiVal;
      }
    }

    // Update previous iteration
    for (int i=0; i<m; i++) {
      uApproxOld[i] = uApprox[i];
    }
  }
  delete uApproxOld;
}

// Solves problem until error is less than a tolerance
void Elliptic::SolveWithTol(const double tol) {

  double *uApproxOld;
  uApproxOld = new double[m];
  for(int i=0; i<m; i++) {
    uApproxOld[i] = (*mFunction).init(mNodes[i]);
    uApprox[i] = (*mFunction).init(mNodes[i]);
  }

  double uVal, psi, psiVal;
  double error = 10;
  int k=0;

  // SOR method stops when the grid norm error is less than a specific tolerance
  // (or stops if iterations become too high)
  while((error >= tol)||(k>=100000)) {

    // Implements SOR method
    for(int i=0; i<m; i++) {

      if (i==0) {
        uVal = (mFvec[0]-(mUpper[0]*uApproxOld[1]))/mDiag[0];
      } else if (i==m-1) {
        uVal = (mFvec[m-1]-(mLower[m-2]*uApprox[m-2]))/mDiag[m-1];
      } else {
        uVal = (mFvec[i]-(mLower[i-1]*uApprox[i-1])-(mUpper[i]*uApproxOld[i+1]))/mDiag[i];
      }
      psi = (*mFunction).psi(mNodes[i]);

      psiVal = (w*uVal)+((1.0-w)*uApproxOld[i]);

      // Updates u values
      if ( psi <psiVal) {
        uApprox[i] = psi;
      } else {
        uApprox[i] = psiVal;
      }
    }

    // Update previous iteration
    for (int i=0; i<m; i++) {
      uApproxOld[i] = uApprox[i];
    }

    // Re-calculates error norm
    double sum = 0;
    for(int i=0; i<n-1; i++) {
      sum = sum + fabs(uApprox[i]-uExact[i]);
    }
    error = sqrt(sum *h);
    k = k+1;

  }
  if (error < tol) {
    std::cout << "\nProcess finished after " << k << " iterations.\n";
  } else {
    std::cout << "\nApproximation did not converge.";
  }

  delete uApproxOld;
}

// Solves problem until approximation converges
void Elliptic::SolveConvergence(const double tol) {
  double *uApproxOld;
  uApproxOld = new double[m];
  for(int i=0; i<m; i++) {
    uApproxOld[i] = (*mFunction).init(mNodes[i]);
    uApprox[i] = (*mFunction).init(mNodes[i]);
  }

  double uVal, psi, psiVal;
  double uDiff = 10;
  int k=0;

  // SOR method stops when the difference in iterations of SOR method is below
  // a tolerance (current u - previous u)
  while((uDiff >= tol)||(k>=100000)) {

    // Implements SOR method
    for(int i=0; i<m; i++) {

      if (i==0) {
        uVal = (mFvec[0]-(mUpper[0]*uApproxOld[1]))/mDiag[0];
      } else if (i==m-1) {
        uVal = (mFvec[m-1]-(mLower[m-2]*uApprox[m-2]))/mDiag[m-1];
      } else {
        uVal = (mFvec[i]-(mLower[i-1]*uApprox[i-1])-(mUpper[i]*uApproxOld[i+1]))/mDiag[i];
      }
      psi = (*mFunction).psi(mNodes[i]);

      psiVal = (w*uVal)+((1.0-w)*uApproxOld[i]);

      // Updates u values
      if ( psi <psiVal) {
        uApprox[i] = psi;
      } else {
        uApprox[i] = psiVal;
      }
    }

    // Finds the difference in iterations for all u values and calculates the
    // grid norm for this
    double sum = 0;
    for(int i=0; i<n-1; i++) {
      sum = sum + fabs(uApprox[i]-uApproxOld[i]);
    }
    uDiff = sqrt(sum *h);
    k = k+1;

    // Update previous iteration
    for (int i=0; i<m; i++) {
      uApproxOld[i] = uApprox[i];
    }

  }
  if (uDiff < tol) {
    std::cout << "\nProcess finished after " << k << " iterations.\n";
  } else {
    std::cout << "\nApproximation did not converge.";
  }

  delete uApproxOld;

}

// Calculates unconstrained solution for the corresponding Q1 problem
void Elliptic::UnconstrainedSol() {
  double* uArray;
  uArray = new double[n-1];

  // Creates dummy vectors
  double *delta, *Gvec;
  delta = new double[n-1];
  Gvec = new double[n-1];
  for(int i=0; i<=n-2; i++) {
  delta[i] = mDiag[i];
  Gvec[i] = mFvec[i];
  }

  // Elimination stage
  for(int i=1; i<=n-2; i++) {
    delta[i] = delta[i] - mUpper[i-1]*(mLower[i-1]/delta[i-1]);
    Gvec[i] = Gvec[i] - Gvec[i-1]*(mLower[i-1]/delta[i-1]);
  }

  // Backsolve
  uArray[n-2] = Gvec[n-2]/delta[n-2];
  for(int i=n-3; i>=0; i--) {
    uArray[i] = ( Gvec[i] - mUpper[i]*uArray[i+1] )/delta[i];
  }

  // Updates unconstrained solution
  for(int i=0; i<n-1; i++){
    uUnconstrained[i] = uArray[i];
  }

  // Deallocates storage
  delete delta;
  delete Gvec;
  delete uArray;
}

// Finds exact solution to the elliptical inequality problem
void Elliptic::FindUExact() {

  for(int i=0; i<n-1; i++){
    // Free boundaries
    double x1 = sqrt(3)/double(5);
    double x2 = 1.0-(sqrt(3)/double(5));

    // Defines functions
    double a = -double(25)/double(3);
    double b1 = double(10)/double(sqrt(3));
    double c1 = 0;
    double b2 = (50.0-(10*sqrt(3)))/double(3);
    double c2 = (-25.0+(10*sqrt(3)))/double(3);

    // u(x) is a piecewise function (depending on the free boundaries)
    if (mNodes[i] < x1) {
      uExact[i] = (a*pow(mNodes[i],2))+(b1*mNodes[i])+c1;
    } else if (mNodes[i] > x2) {
      uExact[i] = (a*pow(mNodes[i],2))+(b2*mNodes[i])+c2;
    } else {
      uExact[i] = (*mFunction).psi(mNodes[i]);
    }
  }
}

// Show approximation
void Elliptic::ShowApprox() {
  std::cout.precision(6);
  std::cout << "\nApproximation: ";
  for (int i=0; i<n-1; i++) {
    std::cout << uApprox[i] << " ";
  }
  std::cout << std::endl;
}

// Show exact solution
void Elliptic::ShowExact() {
  std::cout.precision(6);
  std::cout << "\nExact solution: ";
  for (int i=0; i<n-1; i++) {
    std::cout << uExact[i] << " ";
  }
  std::cout << std::endl;
}

// Shows the grid norm
void Elliptic::ShowNorm() {
  std::cout.precision(6);
  double sum = 0;
  for(int i=0; i<n-1; i++) {
    sum = sum +  fabs(uApprox[i]-uExact[i]);
  }
  sum = sqrt(sum *h);
  std::cout << "\nGrid norm: " << sum << "\n";
}

// Retrieves grid norm
double Elliptic::GetNorm() {
  std::cout.precision(6);
  double sum = 0;
  for(int i=0; i<n-1; i++) {
    sum = sum +  fabs(uApprox[i]-uExact[i]);
  }
  sum = sqrt(sum *h);
  return sum;

}
void Elliptic::PlotApproximation(std::string constraint) {
  // Writes data to a file
  system("rm EllipticIneqPlot.csv");
  std::ofstream file;
  file.open("EllipticIneqPlot.csv");
  assert(file.is_open());

  // Saves x values
  file << 0 << ",";
  for(int i=0; i<n-1; i++) {
    file << mNodes[i] << ",";
  }
  file << 1 << ",";
  file << std::endl;

  // Saves approximation
  file << alpha << ",";
  for(int i=0; i<n-1; i++) {
    file << uApprox[i] << ",";
  }
  file << beta << ",";
  file << std::endl;

  // Saves exact solution or unconstrained solution (depending on 'constraint')
  file << alpha << ",";
  for(int i=0; i<n-1; i++) {
    if (constraint == "unconstrained") {
      file << uUnconstrained[i] << ",";
    } else { // if constraint is "constrained"
      file << uExact[i] << ",";
    }
  }
  file << beta << ",";

  file.close();
  system("cp EllipticIneqPlot.csv ../../../MATLAB/");
}

// Deallocates storage
Elliptic::~Elliptic() {
  delete mNodes;
  delete mDiag;
  delete mUpper;
  delete mLower;
  delete mFvec;
  delete uApprox;
  delete uExact;
  delete uUnconstrained;
}
