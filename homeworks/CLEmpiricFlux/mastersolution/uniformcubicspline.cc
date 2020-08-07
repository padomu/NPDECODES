/**
 * @file uniformcubicspline.cc
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 18.07.2019
 * @copyright Developed at ETH Zurich
 */

#include "uniformcubicspline.h"

#include <Eigen/Core>
#include <cassert>

namespace {

template <typename T>
constexpr T Square(T x) {
  return x * x;
}

template <typename T>
constexpr T Cube(T x) {
  return x * x * x;
}

constexpr int getJ(double a, double b, unsigned int n, double u) {
  return u < b ? (int)(n * ((u - a) / (b - a)) + 1.0) : n;
}

constexpr double zeta(double a, double b, unsigned int n, double j) {
  return a + j * (b - a) / n;
}

}  // namespace

namespace CLEmpiricFlux {

UniformCubicSpline::UniformCubicSpline(double a, double b, Eigen::VectorXd f,
                                       Eigen::VectorXd M)
    : _n(f.size() - 1), _a(a), _b(b), _f(std::move(f)), _M(std::move(M)) {
  assert(b >= a);
  assert(_f.size() >= 2);
  assert(_f.size() == _M.size());
}

double UniformCubicSpline::operator()(double u) const {
  assert((u >= _a) && (u <= _b));
  int j = getJ(_a, _b, _n, u);
  double h = (_b - _a) / _n;
  double tau = (u - zeta(_a, _b, _n, j - 1)) / h;

  return _f(j) * tau + _f(j - 1) * (1.0 - tau) +
         Square(h) / 6.0 *
             (_M(j) * (Cube(tau) - tau) +
              _M(j - 1) * (-Cube(tau) + 3.0 * Square(tau) - 2.0 * tau));
}

double UniformCubicSpline::derivative(double u) const {
  assert((u >= _a) && (u <= _b));
  int j = getJ(_a, _b, _n, u);
  double h = (_b - _a) / _n;
  double tau = (u - zeta(_a, _b, _n, j - 1)) / h;

  return (_f(j) - _f(j - 1) +
          Square(h) / 6.0 *
              (_M(j) * (3.0 * Square(tau) - 1.0) +
               _M(j - 1) * (-3.0 * Square(tau) + 6.0 * tau - 2.0))) *
         (1.0 / h);
}

}  // namespace CLEmpiricFlux
