#ifndef __YIELDSPLINE_H_INCLUDED__
#define __YIELDSPLINE_H_INCLUDED__
#include <cmath>
#include "Spline.h"
#include <iostream>
#include "BondUtilities.h"
#include "CurveFeatures.h"
#include "Date.h"
#include <sstream>//stringstream
#include <unordered_map>

typedef std::vector<SpotValue> YieldCurve;
typedef std::vector<SpotValue> FutureCurve;
typedef std::vector<SpotValue> LiborCurve;
typedef std::vector<SpotValue> SwapCurve;
class YieldSpline {
private:
  YieldCurve yield;
  FutureCurve forwardCurve;
  std::vector<double> theta;
  double maxMaturity;
   std::vector<double> splineX;
    std::vector<double> splineY;
    std::vector<double> splineZ;
    
  double r0;
  //char type[];
public:
    YieldSpline(YieldCurve&, Date&, double);
    YieldSpline();
    
    
   
    
    void getForwardCurve(auto&);
    void getSpotCurve(auto&);
    void computeSimpleSwapSpline(LiborCurve&, SwapCurve&, Date&);
    void computeSimpleSwapSpline(LiborCurve&, SwapCurve&, Date&, double);
    void computeYieldFunction(Date&);
    void computeSimpleFutureSpline(double, FutureCurve&, Date&);
    double Yield(double);//Cubic Spline
    double Forward(double);//Cubic Spline
    double getShortRate();

};
#include "YieldSpline.hpp"

#endif
