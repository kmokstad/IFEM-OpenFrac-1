// $Id$
//==============================================================================
//!
//! \file HermiteInterpolator.C
//!
//! \date Oct 18 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Piece-wise cubic hermite interpolation utilities.
//!
//==============================================================================

#include "HermiteInterpolator.h"
#include "IFEM.h"
#include <fstream>


double HermiteInterpolator::evaluate (double x, int derOrder) const
{
  size_t i = 0;
  while (i+2 < grid.size() && x > grid[i+1])
    ++i;

  double h  = grid[i+1] - grid[i];
  double df = (values[i+1] - values[i]) / h;

  double c2 = -(2.0*derivs[i] - 3.0*df + derivs[i+1]) / h;
  double c3 =  (    derivs[i] - 2.0*df + derivs[i+1]) / (h*h);

  x -= grid[i];
  switch (derOrder) {
  case 0: return values[i] + x*(derivs[i] + x*(c2 + x*c3));
  case 1: return derivs[i] + x*(           2.0*c2 + x*3.0*c3);
  case 2: return                           2.0*c2 + x*6.0*c3;
  case 3: return                                      6.0*c3;
  }

  return 0.0;
}


bool HermiteInterpolator::findMinimum (double& x) const
{
  const double absTol = 1.0e-20;
  const double relTol = 1.0e-12;
  const double epsZero = 1.0e-10;
  const int maxIts = 20;

#ifdef SP_DEBUG
  std::cout <<"\nHermiteInterpolator::findMinimum:";
#endif

  std::ofstream of("val.asc");
  for (size_t i = 0; i < 100; ++i)
    of << this->evaluate(grid.front() + (grid.back()-grid.front())/99*i) << " ";
  of.close();
  std::ofstream of2("der.asc");
  for (size_t i = 0; i < 100; ++i)
    of2 << this->evaluateDeriv(grid.front() + (grid.back()-grid.front())/99*i) <<  " ";
  of2.close();

  // Newton loop to find zeros
  std::vector<double> extrema;
  for (size_t pidx = 0; pidx < grid.size(); ++pidx)
  {
    x = grid[pidx];
    double dx = 1.0;
    size_t its = 0;
    std::cout << "newton loop for grid point " << pidx <<":";
    while ((dx/x > relTol || dx/x < -relTol) &&
           (dx > absTol || dx < -absTol) && its < 100)
    {
      double I  = this->evaluateDeriv(x);
      double I2 = this->evaluateDeriv2(x);
      if (I2 < epsZero && I2 > -epsZero)
        return false; // breakdown - probably a linear function
      dx = I / I2;
      x -= dx;
      std::cout << "\n\t" << dx << std::endl;
      ++its;
    }
    if (its >= maxIts)
      continue;

#ifdef SP_DEBUG
    std::cout <<" f("<< x <<") = "<< this->evaluate(x);
#endif
    if (x >= grid.front() && x <= grid.back() && this->evaluateDeriv2(x) > 0.0)
      extrema.push_back(x);
  }

  if (!extrema.empty()) {
    double val = 1e100;
    for (double& it : extrema) {
      double val2 = this->evaluate(it);
      if (val2 < val) {
        x = it;
        val = val2;
      }
    }

#ifndef SP_DEBUG
    IFEM::cout <<" alpha^* = "<< x << std::endl;
#endif

    return true;
  }

  IFEM::cout <<"\n\t*** No energy minimum found in [u";
  if (grid.front() == 0.0)
    IFEM::cout <<",u";
  else if (grid.front() == -1.0)
    IFEM::cout <<"-du,u";
  else
    IFEM::cout << grid.front() <<"*du,u";
  if (grid.back() == 1.0)
    IFEM::cout <<"+du]";
  else
    IFEM::cout << grid.front() <<"*du]";
  IFEM::cout << std::endl;
  return false;
}
