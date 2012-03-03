// Written by Peter Kutz.

#include "constants.h"

// TODO: Adjust the various settings below, and add more if you want.

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
//const int theDim[3] = {2, 2, 1};
const int theDim[3] = {20, 20, 1};
#else
//const int theDim[3] = {12, 12, 4};
//const int theDim[3] = {2, 2, 1};
//const int theDim[3] = {12, 12, 1};
//const int theDim[3] = {6, 6, 1};
//const int theDim[3] = {4, 4, 1};
//const int theDim[3] = {14, 14, 1};
const int theDim[3] = {40,40,1};
//const int theDim[3] = {40,40,20};
#endif

//const double theCellSize = 0.5;
const double theCellSize = 1.0;

const double fluidDensity = 1.0;

// alpha - density coef
// beta - temperature coef
//const double buoyAlpha = 0.0;
//const double buoyBeta = 0.0;
//const double Tamb = 300.0;
const double buoyAlpha = 0.0002;
const double buoyBeta = 0.010;
const double Tamb = 0.0;
const double Tmin = -100.0;
const double Tmax = 100.0;


//const double vorticityEpsilon = 0.30;
const double vorticityEpsilon = 0.20;
//const double vorticityEpsilon = 0.00;

