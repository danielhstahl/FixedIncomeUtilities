#ifndef __SPLINE_H_INCLUDED__
#define __SPLINE_H_INCLUDED__
#include <cmath>
#include <vector>
/*Given arrays x[0..n] and y[0..n] containing a tabulated function, i.e., yi = f(xi), with
x1 < x2 <...< xN , and given values yp1 and ypn for the first derivative of the interpolating
function at points 0 and n, respectively, this routine returns an array y2[0..n] that contains
the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or
ypn are equal to 1 × 1030 or larger, the routine is signaled to set the corresponding boundary
condition for a natural spline, with zero second derivative on that boundary.*/
std::vector<double> spline(std::vector<double>&, std::vector<double>&);
double splint(std::vector<double>&, std::vector<double>&, std::vector<double>&, double);
double splintD(std::vector<double>&, std::vector<double>&, std::vector<double>&, double);
#include "Spline.hpp"

#endif
