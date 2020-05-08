#include "Function1.hpp"
#include<cmath>

// Constructor
Function1::Function1() {
}

// Specifies f(x) and evaluates it at point x
double Function1::evaluateF(double x) {
  return (-pow(M_PI,2)*sin(M_PI*x));
}

// Specifies the exact u(x)
double Function1::exactU(double x) {
  return (x-sin(M_PI*x));
}

// Destructor
Function1::~Function1() {

}
