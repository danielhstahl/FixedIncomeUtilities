#include "YieldSpline.h"

YieldSpline::YieldSpline(YieldCurve& yield_, Date& currDate, double r0_){
    yield=yield_;
    r0=r0_;
    computeYieldFunction(currDate);

}
YieldSpline::YieldSpline(){
}
double YieldSpline::getShortRate(){
    return r0;
}
void YieldSpline::getForwardCurve(){
    std::cout<<"{\"Forward\":["; //to send to node
    std::cout<<"{\"x\":"<<splineX[1]<<", \"y\":"<<Forward(splineX[1])<<"}";
    int n=splineX.size();
    for(int i=2;i<n;++i){
        std::cout<<", {\"x\":"<<splineX[i]<<", \"y\":"<<Forward(splineX[i])<<"}";
    }
    std::cout<<"]}\\n"<<std::endl;
}
void YieldSpline::getSpotCurve(){
    std::cout<<"{\"Spot\":["; //to send to node
    std::cout<<"{\"x\":"<<splineX[1]<<", \"y\":"<<splineY[1]/splineX[1]<<"}";
    int n=splineX.size();
    for(int i=2;i<n;++i){
        std::cout<<", {\"x\":"<<splineX[i]<<", \"y\":"<<splineY[i]/splineX[i]<<"}";
    }
    std::cout<<"]}\\n"<<std::endl;
}
void YieldSpline::computeYieldFunction(Date& currDate){//Cubic Spline
    int n=yield.size();
  //Date currDate;
    splineX=std::vector<double>(n+1);
    splineY=std::vector<double>(n+1);
    double dt=0;
    splineX[0]=0;
    splineY[0]=0;
    yield[0].date.setScale("year");
    dt=yield[0].date-currDate;
    splineX[1]=dt;
    splineY[1]=convertLiborToContinuous(yield[0].value, dt)*dt; //convert to continuous compounding ASSUMING that yield[i].value is a simple rate (not  compounded at all), that is, 1/(1+yield*t)=bond.
    for(int i=1;i<n;++i){
        yield[i].date.setScale("year");
        dt=yield[i].date-currDate;
        splineX[i+1]=dt;
        splineY[i+1]=convertLiborToContinuous(yield[i].value, dt)*dt; //convert to continuous compounding ASSUMING that yield[i].value is a simple rate (not compounded at all), that is, 1/(1+yield*t)=bond.
  }
  splineZ=spline(splineX, splineY, r0);//, .03777586); //TEMPOORARY
  /*std::cout<<"{\"Spot\":";
  getSpotCurve();
  std::cout<<"}"<<std::endl;
  std::cout<<"{\"Forward\":";
  getForwardCurve();
  std::cout<<"}"<<std::endl;*/
}
void YieldSpline::computeSimpleSwapSpline(LiborCurve &libor, SwapCurve &swap, Date& currDate){
    computeSimpleSwapSpline(libor, swap, currDate, .25);
}
void YieldSpline::computeSimpleSwapSpline(LiborCurve &libor, SwapCurve &swap, Date& currDate, double delta){
  currDate.setScale("year");
    int liborSize=libor.size();
    libor[0].date.setScale("year");
    r0=convertLiborToContinuous(libor[0].value, libor[0].date-currDate);
  //get libor for every three month interval (0, 3, 9, 12)
    std::vector<double> splineLiborX(liborSize);
    std::vector<double> splineLiborY(liborSize);
    double dt=0;
    std::unordered_map<int, double> trackLibor;
    int integerTime=0;
    trackLibor.insert({0, 1.0});
    for(int i=0; i<liborSize; ++i){
        libor[i].date.setScale("year");
        dt=libor[i].date-currDate;
        yield.push_back(SpotValue(libor[i].date, convertLiborToContinuous(libor[i].value, dt)));
        integerTime=(int)(dt/delta);
    //if(trackLibor.find(integerTime)==trackLibor.end()){
        trackLibor.insert({integerTime, convertLiborToBond(libor[i].value, dt)});
    //}

        splineLiborX[i]=dt;
        splineLiborY[i]=libor[i].value;
    }
    double getSumZeroBonds=0;
    std::vector<double> splineLiborZ=spline(splineLiborX, splineLiborY);
    for(int i=1; i <=integerTime; ++i){
        double val=0;
        dt=(double)i*delta;
    //std::cout<<"Time: "<<dt;
        if(trackLibor.find(i)==trackLibor.end()){
            val=splint(splineLiborX, splineLiborY, splineLiborZ, dt);//missing
            val=convertLiborToBond(val, dt);
      //std::cout<<" and we are adding a libor! its value is "<<val<<std::endl;
        }
        else{
            val=trackLibor.find(i)->second;
      //std::cout<<" and we already had a libor and its value is "<<val<<std::endl;
        }
        getSumZeroBonds+=val;
    }

  //get swap curve
    int swapSize=swap.size();
    std::vector<double> splineSwapX(swapSize);
    std::vector<double> splineSwapY(swapSize);
    dt=0;
    for(int i=0; i<swapSize; ++i){
        swap[i].date.setScale("year");
        dt=swap[i].date-currDate;
        splineSwapX[i]=dt;
        splineSwapY[i]=swap[i].value;
    }
    std::vector<double> splineSwapZ=spline(splineSwapX, splineSwapY);
    int maxTime=(int)dt/delta;
  //get implied libor curve
    for(int i=integerTime+1; i<maxTime;i++){
        dt=delta*i;
        double swpYield=splint(splineSwapX, splineSwapY, splineSwapZ, delta*(i-1));
        double impliedBond=(1.0-swpYield*delta*getSumZeroBonds)/(swpYield*delta+1.0);
        getSumZeroBonds+=impliedBond;
        yield.push_back(SpotValue(currDate+dt, convertLiborToContinuous(convertBondToLibor(impliedBond, dt), dt)));
    }
    computeYieldFunction(currDate);
}
void YieldSpline::computeSimpleFutureSpline(double r0_, FutureCurve &future, Date& currDate){ //converts future curve into continously compounded zero coupon curve.
  int n=future.size();
  r0=r0_;
  //Date currDate;
  std::vector<double> splineXF(n+1);
  std::vector<double> splineYF(n+1);
  double dt=0;
  splineXF[0]=0;
  splineYF[0]=0;
  //std::cout<<n<<std::endl;
  for(int i=0; i<n; ++i){
    future[i].date.setScale("year");
    dt=future[i].date-currDate;
    splineXF[i+1]=dt;
    splineYF[i+1]=1-future[i].value;//-vs.EuroDollarFutureAdjustment(r0, 0, .25, dt)*.25;
  }
  int m=(int)dt*4;
  std::vector<double> splineZF=spline(splineXF, splineYF, r0);
  double bond=1;
  for(int i=0; i<m; ++i){
    dt=.25*(i+1);
    bond=bond/(.25*splint(splineXF, splineYF, splineZF, dt)+1);
    yield.push_back(SpotValue(currDate+dt, (1-bond)/dt));
  }
  computeYieldFunction(currDate);
}
double YieldSpline::Yield(double t){//Cubic Spline
  return splint(splineX, splineY, splineZ, t);
}

double YieldSpline::Forward(double t){//Cubic Spline
  return splintD(splineX, splineY, splineZ, t);
}
