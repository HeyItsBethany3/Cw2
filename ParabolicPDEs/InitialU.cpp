#include "InitialU.hpp"
#include<cmath>

// Constructor
InitialU::InitialU() {
}

// Specifies initial u_0(x) and evaluates it at point x
double InitialU::evaluate(double x) {
  return sin(2*M_PI*x)+x;
}


// Destructor
InitialU::~InitialU() {

}
