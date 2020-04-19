#include "Functions.hpp"
#include<cmath>

// Constructor
Functions::Functions(const double strike, const double interest, const double sigma) {
  K = strike;
  r = interest;
  vol = sigma;
}

double Functions::payoff(double x) {
  double payoff;
  if (K>x) {
    payoff = K-x;
  } else {
    payoff = 0;
  }
  return payoff;

}

// Destructor
Functions::~Functions() {

}
