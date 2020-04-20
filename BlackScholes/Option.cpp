#include "Option.hpp"
#include <cmath>
#include <iostream>

// Constructor
Option::Option(const double strike, const double interest, const double sigma,
      const double maturity, const double maxX,
      AbstractFunctions& aFunction, const int N, const int L) {

        K = strike;
        r = interest;
        vol = sigma;
        T = maturity;
        R = maxX;
        mFunction = &aFunction;

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
        uApprox = new double[m];

}

// Destructor
Option::~Option() {
  delete xNodes;
  delete mDiag;
  delete mUpper;
  delete mLower;
  delete uApprox;
}

// Construct matrix A
void Option::ConstructMatrix() {


  // Diagonal elements of A
  for (int i=0; i<m; i++) {
    double lamda = xNodes[i]/h;
    mDiag[i] = 1+(r*deltaT)+(r*deltaT*lamda)+(pow(vol,2)*pow(lamda,2));
  }

  // Construct upper diagonal elements of A
  for (int i=0; i<m-1; i++) {
    double lamda = xNodes[i]/h;
    mUpper[i]=-((pow(vol,2)*pow(lamda,2))/double(2))-(r*lamda);
  }

  // Construct lower diagonal elements of A
  for (int i=1; i<m; i++) {
    double lamda = xNodes[i]/h;
    mLower[i]=-(pow(vol,2)*pow(lamda,2)*deltaT)/double(2);
  }
}

// Displays system to solve
void Option::ShowMatrix() {

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

// Solve problem
void Option::Approximate() {
  double t = deltaT; // Is this right?
  double* uApproxOld;
  uApproxOld = new double[m];
  double* uApproxNew;
  uApproxNew = new double[m];

  for(int j=0; j<m; j++) {
    uApproxOld[j] = (*mFunction).payoff(xNodes[j]);
  }
  double factor1 = -((pow(vol,2)*pow((xNodes[0]/h),2)*deltaT)/double(2));
  uApproxOld[0] += -(factor1 * (*mFunction).f0(t));
  double factor2 = -((pow(vol,2)*pow((xNodes[m-1]/h),2))/double(2))-(r*(xNodes[m-1]/h));
  uApproxOld[m-1] += -(factor2 * (*mFunction).exactU(xNodes[m-1], t));

  for(int i=1; i<=l; i++) {

    // Solve tridiagonal system of equations A u_n+1 = u_n
    double *delta;
    delta = new double[n-1];

    for(int i=0; i<=n-2; i++) {
      delta[i] = mDiag[i];
    }

    // Elimination stage
    for(int i=1; i<=n-2; i++)
    {
      delta[i] = delta[i] - mUpper[i-1]*(mLower[i-1]/delta[i-1]);
      uApproxOld[i] = uApproxOld[i] - uApproxOld[i-1]*(mLower[i-1]/delta[i-1]);
    }

    // Backsolve
    uApproxNew[n-2] = uApproxOld[n-2]/delta[n-2];
    for(int i=n-3; i>=0; i--)
    {
      uApproxNew[i] = ( uApproxOld[i] - mUpper[i]*uApproxNew[i+1] )/delta[i];
    }

    // Deallocates storage
    delete delta;

    // Update time
    t += deltaT;

    // Update old vector for next iteration (and add boundary conditions)
    for(int i=0; i<m; i++) {
      uApproxOld[i] = uApproxNew[i];
    }
    double factor1 = -((pow(vol,2)*pow((xNodes[0]/h),2)*deltaT)/double(2));
    uApproxOld[0] += -(factor1 * (*mFunction).f0(t));
    double factor2 = -((pow(vol,2)*pow((xNodes[m-1]/h),2))/double(2))-(r*(xNodes[m-1]/h));
    uApproxOld[m-1] += -(factor2 * (*mFunction).exactU(xNodes[m-1], t));

  }

  // Save uApprox
  for(int i=0; i<m; i++) {
    uApprox[i] = uApproxNew[i];
  }


  delete uApproxOld;
  delete uApproxNew;
}

// Show approximation
void Option::ShowApprox() {
  std::cout << "Approx: ";
  for(int i=0; i<m; i++) {
    std::cout << uApprox[i] << " ";
  }
  std::cout << std::endl;
}

void Option::ShowExact() {
  std::cout << "Exact: ";
  for(int i=0; i<m; i++) {
    std::cout << (*mFunction).exactU(xNodes[i],T) << " ";
  }
  std::cout << std::endl;
}

void Option::ShowError() {
  std::cout << "Error: ";
  for(int i=0; i<m; i++) {
    std::cout << fabs(uApprox[i]-(*mFunction).exactU(xNodes[i], T)) << " ";
  }
  std::cout << std::endl;
}

// Shows the grid norm
void Option::Norm() {
  double sum = 0;
  for(int i=0; i<n-1; i++) {
    sum = sum +  fabs(uApprox[i]-(*mFunction).exactU(xNodes[i], T));
  }
  sum = sqrt(sum *h);
  std::cout << "\nGrid norm: " << sum << "\n";

}
