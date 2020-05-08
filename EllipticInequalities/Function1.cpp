#include "Function1.hpp"
#include<cmath>

// Constructor
Function1::Function1() {
}

// Specifies f(x) and evaluates it at point x
double Function1::f(double x) {
  return double(50)/double(3);
}

// psi
double Function1::psi(double x) {
  return 1;
}

// Initial guess for u at all interior points
double Function1::init(double x) {
  return 0;
}

// Destructor
Function1::~Function1() {

}
