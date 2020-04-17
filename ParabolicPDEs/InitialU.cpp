#include "InitialU.hpp"
#include<cmath>

// Constructor
InitialU::InitialU() {
}

// Specifies f(x) and evaluates it at point x
double InitialU::evaluate(double x) {
  return sin(2*M_PI*x)+x;
}


// Destructor
InitialU::~InitialU() {

}
