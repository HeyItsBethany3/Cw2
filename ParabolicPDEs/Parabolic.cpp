#include "Parabolic.hpp"
#include <iostream>
#include <cmath>
#include <fstream>


// Constructor
Parabolic::Parabolic(const double A, const double time, const double g_0,
  const double g_1, Function1D& InitialU, Function2D& ExactU,
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
    // Calculates nodes
    for(int i=0; i<m; i++) {
      xNodes[i]=(i+1)*h;
    }

    mDiag = new double[m];
    mUpper = new double[m-1];
    mLower = new double[m-1];
    uApprox = new double[m];
}

// Construct matrix A
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

// Displays matrix A
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

// Finds approximation for u(T,x)
void Parabolic::Approximate() {

  double* uApproxOld;
  uApproxOld = new double[m];
  double* uApproxNew;
  uApproxNew = new double[m];
  double lamda = ((a*deltaT)/double(pow(h,2)));

  // Use initial condition u(x,0)
  for(int j=0; j<m; j++) {
    uApproxOld[j] = (*mInitialU).evaluate(xNodes[j]);
  }
  // Add boundary information
  uApproxOld[0] += lamda*g0;
  uApproxOld[m-1] += lamda*g1;

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

    delete delta;

    // Update old vector for next iteration (and add boundary information)
    for(int i=0; i<m; i++) {
      uApproxOld[i] = uApproxNew[i];
    }
    uApproxOld[0] += lamda*g0;
    uApproxOld[m-1] += lamda*g1;

  }
  // Save approximation for u(x,T)
  for(int i=0; i<m; i++) {
    uApprox[i] = uApproxNew[i];
  }
  delete uApproxOld;
  delete uApproxNew;
}

// Show approximation
void Parabolic::ShowApprox() {
  std::cout << "Approx: ";
  for(int i=0; i<m; i++) {
    std::cout << uApprox[i] << " ";
  }
  std::cout << std::endl;
}

// Shows exact value at T
void Parabolic::ShowExact() {
  std::cout << "Exact: ";
  for(int i=0; i<m; i++) {
    std::cout << (*mExactU).evaluate(xNodes[i],T) << " ";
  }
  std::cout << std::endl;
}

 // Show absolute errors
void Parabolic::ShowError() {
  std::cout << "Error: ";
  for(int i=0; i<m; i++) {
    std::cout << fabs(uApprox[i]-(*mExactU).evaluate(xNodes[i], T)) << " ";
  }
  std::cout << std::endl;
}

// Shows the max error norm
void Parabolic::ShowNorm() {
  double sum = 0;
  for(int i=0; i<n-1; i++) {
    sum = sum +  fabs(uApprox[i]-(*mExactU).evaluate(xNodes[i], T));
  }
  sum = sqrt(sum *h);
  std::cout << "\nError norm: " << sum << "\n";
}

// Retrieves max error norm
double Parabolic::GetMaxError() {
  double error = 0;
  for(int i=0; i<n-1; i++) {
    double ei = fabs(uApprox[i]-(*mExactU).evaluate(xNodes[i], T));
    if (ei > error) {
      error = ei;
    }
  }
  return error;
}

// Saves nodes and initial u(x,0)
void Parabolic::SaveInitial() {
  std::ofstream file;
  file.open("ParabolicPlot.csv", std::ios::app);
  assert(file.is_open());

  // x values
  file << 0 << ",";
  for(int i=0; i<n-1; i++) {
    file << xNodes[i] << ",";
  }
  file << 1 << ",";
  file << std::endl;

  // u at time 0
  file << (*mInitialU).evaluate(0) << ",";
  for(int i=0; i<n-1; i++) {
    file << (*mInitialU).evaluate(xNodes[i]) << ",";
  }
  file << (*mInitialU).evaluate(1) << ",";
  file << std::endl;
  file.close();
}

// Saves approximation and exact solution at T
void Parabolic::SaveApprox() {
  std::ofstream file;
  file.open("ParabolicPlot.csv", std::ios::app);
  assert(file.is_open());

  // Approximation
  file << g0 << ",";
  for(int i=0; i<n-1; i++) {
    file << uApprox[i] << ",";
  }
  file << g1 << ",";
  file << std::endl;

  // Exact solution
  file << g0 << ",";
  for(int i=0; i<n-1; i++) {
    file << (*mExactU).evaluate(xNodes[i], T) << ",";
  }
  file << g1 << ",";
  file << std::endl;
  file.close();
}


// Destructor
Parabolic::~Parabolic() {
  delete xNodes;
  delete mDiag;
  delete mUpper;
  delete mLower;
  delete uApprox;
}
