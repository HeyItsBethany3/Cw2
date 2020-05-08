#include "Option.hpp"
#include <cmath>
#include <iostream>
#include <fstream>

// Constructor
Option::Option(const double strike, const double interest, const double sigma,
      const double maturity, const double maxX, const double weight,
      AbstractFunctions& aFunction, const int N, const int L) {

        K = strike;
        r = interest;
        vol = sigma;
        T = maturity;
        R = maxX;
        mFunction = &aFunction;
        w = weight;

        n=N;
        m=N-1;
        l=L;
        deltaT = T/double(l); // time step-size

        xNodes = new double[m];
        h = R/double(N); // spatial step-size
        for(int i=0; i<m; i++) {
          xNodes[i]=(i+1)*h;
        }

        mDiag = new double[m];
        mUpper = new double[m-1];
        mLower = new double[m-1];
        uApprox = new double[m];
        European = new double[m];
        FreeBoundary = new double[l];

}

// Destructor
Option::~Option() {
  delete xNodes;
  delete mDiag;
  delete mUpper;
  delete mLower;
  delete uApprox;
  delete European;
  delete FreeBoundary;
}

// Construct matrix A
void Option::ConstructMatrix() {

  double val;

  // Diagonal elements of A
  for (int i=0; i<m; i++) {
    val = 1+(r*deltaT)+((r*deltaT*xNodes[i])/double(h));
    val = val + ((pow(vol,2)*pow(xNodes[i],2)*deltaT)/double(pow(h,2)));
    mDiag[i] = val;
  }

  // Construct upper diagonal elements of A
  for (int i=0; i<m-1; i++) {
    val = -((pow(vol,2)*pow(xNodes[i],2)*deltaT)/double(pow(h,2)*2));
    val += -((r*xNodes[i]*deltaT)/double(h));
    mUpper[i] = val;
  }

  // Construct lower diagonal elements of A
  for (int i=1; i<m; i++) {
    mLower[i-1]= -((pow(vol,2)*pow(xNodes[i],2)*deltaT)/double(2*pow(h,2)));
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

// European approximation
void Option::FindEuropean() {

  double t = deltaT;
  double* uApproxOld;
  uApproxOld = new double[m];
  double* uApproxNew;
  uApproxNew = new double[m];

  for(int j=0; j<m; j++) {
    uApproxOld[j] = (*mFunction).payoff(xNodes[j]);
  }
  double factor1 = -((pow(vol,2)*pow(xNodes[0],2)*deltaT)/double(2*pow(h,2)));
  uApproxOld[0] += -(factor1 * (*mFunction).f0(t));
  double factor2 = -((pow(vol,2)*pow(xNodes[m-1],2)*deltaT)/double(pow(h,2)*2));
  factor2 += -((r*xNodes[m-1]*deltaT)/double(h));
  uApproxOld[m-1] += -(factor2 * (*mFunction).fR(R, t));

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
    double factor1 = -((pow(vol,2)*pow(xNodes[0],2)*deltaT)/double(2*pow(h,2)));
    uApproxOld[0] += -(factor1 * (*mFunction).f0(t));
    double factor2 = -((pow(vol,2)*pow(xNodes[m-1],2)*deltaT)/double(pow(h,2)*2));
    factor2 += -((r*xNodes[m-1]*deltaT)/double(h));
    uApproxOld[m-1] += -(factor2 * (*mFunction).fR(R, t));

  }

  // Save uApprox
  for(int i=0; i<m; i++) {
    European[i] = uApproxNew[i];
  }

  delete uApproxOld;
  delete uApproxNew;
}


// iter is iterations per time step
void Option::SolveWithIter(const int iter) {

  double *uApproxOld;
  uApproxOld = new double[m];
  for(int i=0; i<m; i++) {
    uApproxOld[i] = (*mFunction).payoff(xNodes[i]);
    uApprox[i] = (*mFunction).payoff(xNodes[i]);
  }

  double uVal, payoff, psiVal;
  double* fArray;
  fArray = new double[m];
  double t = deltaT;
  for(int i=1; i<=l; i++) {

    // Update f Array
    for(int i=0; i<m; i++) {
      fArray[i] = uApproxOld[i];
    }
    double factor1 = -((pow(vol,2)*pow(xNodes[0],2)*deltaT)/double(2*pow(h,2)));
    fArray[0] += -(factor1 * (*mFunction).f0(t));
    double factor2 = -((pow(vol,2)*pow(xNodes[m-1],2)*deltaT)/double(pow(h,2)*2));
    factor2 += -((r*xNodes[m-1]*deltaT)/double(h));
    fArray[m-1] += -(factor2 * (*mFunction).fR(R,t));

    for(int k=1; k<=iter; k++) {
      for(int i=0; i<m; i++) {

        if (i==0) {
          uVal = (fArray[0]-(mUpper[0]*uApproxOld[1]))/mDiag[0];
        } else if (i==m-1) {
          uVal = (fArray[m-1]-(mLower[m-2]*uApprox[m-2]))/mDiag[m-1];
        } else {
          uVal = (fArray[i]-(mLower[i-1]*uApprox[i-1])-(mUpper[i]*uApproxOld[i+1]))/mDiag[i];
        }
        payoff = (*mFunction).payoff(xNodes[i]);
        psiVal = (w*uVal)+((1.0-w)*uApproxOld[i]);

        if ( payoff >= psiVal) {
          uApprox[i] = payoff;
        } else {
          uApprox[i] = psiVal;
        }
      }
      for (int i=0; i<m; i++) {
        uApproxOld[i] = uApprox[i];
      }
    }

    // Find free boundary x value for this time step
    bool FBNotFound = true;
    int k=0;
    double FB = R+1;  // Represents infinity stopping time
    double payoff;
    while((FBNotFound)&&(k<l)) {
      payoff =  (*mFunction).payoff(xNodes[k]);
      if (fabs(uApprox[i]-payoff)<10e-10) {
        FB = xNodes[k];
        FBNotFound = false;
      }
      k++;
    }
    FreeBoundary[i-1] = FB;

    // Update time
    t += deltaT;


    for (int i=0; i<m; i++) {
      fArray[i] = 0;
    }
  }


  delete fArray;
  delete uApproxOld;
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
  for(int i=0; i<5; i++) {
    std::cout << fabs(uApprox[i]-(*mFunction).exactU(xNodes[i], T)) << " ";
  }
  std::cout << std::endl;
}

// Shows the grid norm
void Option::ShowNorm() {
  double sum = 0;
  for(int i=0; i<n-1; i++) {
    sum = sum +  fabs(uApprox[i]-(*mFunction).exactU(xNodes[i], T));
  }
  sum = sqrt(sum *h);
  std::cout << "\nGrid norm: " << sum << "\n";

}

double Option::GetMaxError() {
  double error = 0;
  for(int i=0; i<n-1; i++) {
    double ei = fabs(uApprox[i]-(*mFunction).exactU(xNodes[i],T));
    if (ei > error) {
      error = ei;
    }
  }
  return error;
}

void Option::SaveInitial() {
  std::ofstream file;
  file.open("BSIneqPlot.csv", std::ios::app);
  assert(file.is_open());

  // x values
  for(int i=0; i<n-1; i++) {
    file << xNodes[i] << ",";
  }
  file << std::endl;

  // u at tiime 0
  for(int i=0; i<n-1; i++) {
    file << (*mFunction).payoff(xNodes[i]) << ",";
  }
  file << std::endl;
  file.close();
}

void Option::SaveApprox() {
  std::ofstream file;
  file.open("BSIneqPlot.csv", std::ios::app);
  assert(file.is_open());

  // Approximation
  for(int i=0; i<n-1; i++) {
    file << uApprox[i] << ",";
  }
  file << std::endl;

  // European option price
  for(int i=0; i<n-1; i++) {
    file << European[i] << ",";
  }
  file << std::endl;
  file.close();
}

void Option::SaveFB() {

  std::ofstream file;
  file.open("BSIneqFB.csv", std::ios::app);
  assert(file.is_open());

  // Approximation
  for(int i=1; i<=l; i++) {
    file << (i*deltaT) << ",";
  }
  file << std::endl;

  // European option price
  for(int i=0; i<l; i++) {
    file << FreeBoundary[i] << ",";
  }
  file << std::endl;
  file.close();

}
