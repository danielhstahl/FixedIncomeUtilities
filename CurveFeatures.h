#ifndef __MARKETDATA_H_INCLUDED__
#define __MARKETDATA_H_INCLUDED__
#include "Date.h"
#include <vector>
 struct AssetFeatures{
    Date Maturity;
    Date UnderlyingMaturity;
    double Strike;
    double Tenor;
    std::vector<Date> Coupons;
    double CouponRate;
    int type; //types are defined in HullWhiteEngine
     
 };
 struct ForwardValue{
   Date beginDate;
   Date endDate;
   double value;
   ForwardValue(const Date &dt1, const Date &dt2, double val){
     beginDate=dt1;
     endDate=dt2;
     value=val;
   }
 };

 struct SpotValue{
   Date date;
   double value;
   SpotValue(const Date &dt, double val){
     date=dt;
     value=val;
   }
 };
#endif
