#include "Function2.hpp"
#include<cmath>

// Constructor
Function2::Function2() {
}

// Specifies f(x) and evaluates it at point x
double Function2::evaluate(double x) {
  return (exp(-4*T)*sin(2*M_PI*x))+x;
}

// Change time step T
void Function2::changeT(const double t) {
  T = t;
}

// Destructor
Function2::~Function2() {

}
