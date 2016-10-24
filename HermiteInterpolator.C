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
#include <sstream>
#include <cmath>


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
#ifdef SP_DEBUG
  std::cout <<"\nHermiteInterpolator::findMinimum:";
#endif

  static int lsidx = 1;
  std::stringstream str;
  str << "val" << lsidx << ".asc";
  const size_t samples = 200;
  std::ofstream of(str.str());
  for (size_t i = 0; i < samples; ++i)
    of << this->evaluate(grid.front() + (grid.back()-grid.front())/(samples-1)*i) << " ";
  of.close();
  str.str("");
  str << "der" << lsidx << ".asc";
  std::ofstream of2(str.str());
  for (size_t i = 0; i < samples; ++i)
    of2 << this->evaluateDeriv(grid.front() + (grid.back()-grid.front())/(samples-1)*i) <<  " ";
  of2.close();

  // Loop to find zeros
  std::vector<double> extrema;
  for (size_t pidx = 0; pidx < grid.size()-1; ++pidx) {
    double h  = grid[pidx+1] - grid[pidx];
    double df = (values[pidx+1] - values[pidx]) / h;

    double c2 = -(2.0*derivs[pidx] - 3.0*df + derivs[pidx+1]) / h;
    double c3 =  (    derivs[pidx] - 2.0*df + derivs[pidx+1]) / (h*h);

    double a = 3*c3;
    double b = 2*c2 - 3*c3*2*grid[pidx];
    double c = derivs[pidx] - 2*c2*grid[pidx] + 3*c3*grid[pidx]*grid[pidx];

    // check for complex roots
    double det = b*b - 4*a*c;
    if (det < 0)
      continue;

    double x1 = (-b + sqrt(det))/(2.0*a);
    double x2 = (-b - sqrt(det))/(2.0*a);

#ifdef SP_DEBUG
    std::cout <<" f("<< x1 <<") = "<< this->evaluate(x1);
    std::cout <<" f("<< x2 <<") = "<< this->evaluate(x2);
#endif
    if (x1 >= grid.front() && x1 <= grid.back() && this->evaluateDeriv2(x1) > 0.0)
      extrema.push_back(x1);
    else if (x2 >= grid.front() && x2 <= grid.back() && this->evaluateDeriv2(x2) > 0.0)
      extrema.push_back(x2);
  }

  if (!extrema.empty()) {
    double val = 1e100;
    for (double& it : extrema) {
      double val2 = this->evaluate(it);
      std::cout << "x=" << it << ": " << val2 << std::endl;
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
