// Written by Peter Kutz.

#ifndef CONSTANTS_H
#define CONSTANTS_H



#define LERP(a,b,t) (1-t)*a + t*b

#define EPS (10e-20)

// Don't try to modify the values of these here.
// Modify the values of these in constants.cpp instead.
extern const int theMillisecondsPerFrame;
extern const int theDim[3];
extern const double theCellSize;
extern const double fluidDensity;
extern const double buoyAlpha;
extern const double buoyBeta;
extern const double Tamb; 
extern const double Tmin; 
extern const double Tmax; 
extern const double vorticityEpsilon; 

#endif
