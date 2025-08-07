
// Standard Template Library
#include <complex>
#include <cstdlib>
#include <iostream>
#include "hyp2f1.hpp"
#include "AEAE/complex_functions.H"
#include "AEAE/hyp_2F1.cpp"

using namespace std;

double hyp2f1a(double beta, double z)
{

  return (hyp_2F1(std::complex<double>(1.0,                0.0),
                  std::complex<double>(1.0 / beta,         0.0),
                  std::complex<double>(1.0 + (1.0 / beta), 0.0),
                  std::complex<double>(z,                  0.0) )).real();
}

double hyp2f1b(double beta, double z)
{

  return (hyp_2F1(std::complex<double>(1.0,                0.0),
                  std::complex<double>(2.0 / beta,         0.0),
                  std::complex<double>(1.0 + (2.0 / beta), 0.0),
                  std::complex<double>(z,                  0.0) )).real();
}

double hyp2f1c(double gamma, double z)
{

  return (hyp_2F1(std::complex<double>(1.0,                0.0),
                  std::complex<double>(gamma,         0.0),
                  std::complex<double>(1 + gamma, 0.0),
                  std::complex<double>(z,                  0.0) )).real();
}

double hyp2f1d(double gamma, double z)
{

  return (hyp_2F1(std::complex<double>(1.0,                0.0),
                  std::complex<double>(gamma - 1,          0.0),
                  std::complex<double>(gamma,              0.0),
                  std::complex<double>(z,                  0.0) )).real();
}