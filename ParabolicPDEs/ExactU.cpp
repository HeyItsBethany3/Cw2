#include "ExactU.hpp"
#include<cmath>

// Constructor
ExactU::ExactU() {
}

// Specifies u(x) and evaluates it at point x
double ExactU::evaluate(double x, double T) {
  return (exp(-4*T)*sin(2*M_PI*x))+x;
}

// Destructor
ExactU::~ExactU() {

}
