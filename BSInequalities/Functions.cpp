#include "Functions.hpp"
#include<cmath>
#include<iostream>

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

double Functions::exactU(double x, double t) {
  double d1 = (log(x/K)+((r+(0.5*pow(vol,2)))*t))/double(vol*sqrt(t));
  double d2 = (log(x/K)+((r-(0.5*pow(vol,2)))*t))/double(vol*sqrt(t));
  double price = K*exp(-r*t)*(std::erfc(d2/std::sqrt(2))/double(2));
  price += -(x*(std::erfc(d1/std::sqrt(2))/double(2)));
  return price;
}

double Functions::f0(double t) {
  return K;
}

double Functions::fR(double R, double t) {
  return 0;
}

// Destructor
Functions::~Functions() {

}
