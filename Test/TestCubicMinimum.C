//==============================================================================
//!
//! \file TestCubicMinimum.C
//!
//! \date Jul 13 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for class finding the minimum of a cubic hermite interpolant.
//!
//==============================================================================

#include "HermiteInterpolator.h"

#include "gtest/gtest.h"


class CubicFunction {
public:
  static double value(double x) { return (x*x - x)*x; }
  static double tangent(double x) { return (3.0*x - 2.0)*x; } // 3x^2 - 2x
};


class LinearFunction {
public:
  static double value(double x) { return x; }
  static double tangent(double) { return 1.0; }
};


class QuadraticFunction {
public:
  static double value(double x) { return (x+0.25)*(x+0.75); }
  static double tangent(double x) { return x+x+1.0; }
};


class CubicMinimum {
public:
  static bool Find(double& alpha,
                   const std::vector<double>& params,
                   const std::vector<double>& vals,
                   const std::vector<double>& tgts)
  {
    HermiteInterpolator h(params,vals,tgts);
    return h.findMinimum(alpha);
  }
};


TEST(TestCubicMinimum, CubicFunction)
{
  double alpha;
  std::vector<double> params(10), vals(10), tgts(10);

  for (size_t i = 0; i < 10; ++i) {
    params[i] = -1.0 + 2.0/9.0 * i;
    vals[i] = CubicFunction::value(params[i]);
    tgts[i] = CubicFunction::tangent(params[i]);
  }
  ASSERT_FALSE(CubicMinimum::Find(alpha, params, vals, tgts));

  for (size_t i = 0; i < 10; ++i) {
    params[i] = 1.0/9.0 * i;
    vals[i] = CubicFunction::value(params[i]);
    tgts[i] = CubicFunction::tangent(params[i]);
  }
  ASSERT_TRUE(CubicMinimum::Find(alpha, params, vals, tgts));
  ASSERT_FLOAT_EQ(alpha, 2.0/3.0);
}


TEST(TestCubicMinimum, QuadraticFunction)
{
  double alpha;
  std::vector<double> params(10), vals(10), tgts(10);

  for (size_t i = 0; i < 10; ++i) {
    params[i] = -1.0 + 2.0/9.0 * i;
    vals[i] = QuadraticFunction::value(params[i]);
    tgts[i] = QuadraticFunction::tangent(params[i]);
  }
  ASSERT_TRUE(CubicMinimum::Find(alpha, params, vals, tgts));
  ASSERT_FLOAT_EQ(alpha, -0.5);

  for (size_t i = 0; i < 10; ++i) {
    params[i] = 1.0/9.0 * i;
    vals[i] = QuadraticFunction::value(params[i]);
    tgts[i] = QuadraticFunction::tangent(params[i]);
  }
  ASSERT_FALSE(CubicMinimum::Find(alpha, params, vals, tgts));
}


TEST(TestCubicMinimum, LinearFunction)
{
  double alpha;
  std::vector<double> params(10), vals(10), tgts(10);

  for (size_t i = 0; i < 10; ++i) {
    params[i] = -1.0 + 2.0/9.0 * i;
    vals[i] = LinearFunction::value(params[i]);
    tgts[i] = LinearFunction::tangent(params[i]);
  }
  ASSERT_FALSE(CubicMinimum::Find(alpha, params, vals, tgts));

  for (size_t i = 0; i < 10; ++i) {
    params[i] = 1.0/9.0 * i;
    vals[i] = LinearFunction::value(params[i]);
    tgts[i] = LinearFunction::tangent(params[i]);
  }
  ASSERT_FALSE(CubicMinimum::Find(alpha, params, vals, tgts));
}


TEST(TestCubicMinimum, DiscreteValues)
{
  double alpha;
  std::vector<double> params(21);
  for (size_t i = 0; i < 21; ++i)
    params[i] = -1.0 + 2.0/20.0* i;

  std::vector<double> vals({ 1.62983,1.62567,1.62184,1.61831,1.61509,1.61219,1.60959,1.60731,1.60534,1.60367,1.60232,1.60128,1.60055,1.60013,1.60003,1.60024,1.60076,1.60157,1.60265,1.60403,1.60572});
      
  std::vector<double> tgts({ -0.0564113,-0.0524099,-0.0484112,-0.0444145,-0.0404189,-0.036424,-0.032429,-0.0284333,-0.024436,-0.0204368,-0.0164349,-0.0124297,-0.00842082,-0.00440755,-0.00038872,0.00363645,0.00766877,0.01172,0.015798,0.0198927,0.0240022});

  ASSERT_TRUE(CubicMinimum::Find(alpha, params, vals, tgts));
  ASSERT_NEAR(alpha, -1.0, 1.0e-8);
}
