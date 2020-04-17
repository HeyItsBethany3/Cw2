#include "Function1.hpp"
#include<cmath>

// Constructor
Function1::Function1() {
}

// Specifies f(x) and evaluates it at point x
double Function1::evaluate(double x) {
  return sin(2*M_PI*x)+x;
}


// Destructor
Function1::~Function1() {

}
