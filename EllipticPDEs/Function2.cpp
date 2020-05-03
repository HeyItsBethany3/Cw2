#include "Function2.hpp"
#include<cmath>

// Constructor
Function2::Function2() {
}

// Specifies f(x) and evaluates it at point x
double Function2::evaluateF(double x) {
  return double(50)/double(3);
}

double Function2::exactU(double x) {
  return 0; // No information available
}

// Destructor
Function2::~Function2() {

}
