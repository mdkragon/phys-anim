// Written by Peter Kutz.

#include "constants.h"

// TODO: Adjust the various settings below, and add more if you want.

const int theMillisecondsPerFrame = 10;

#ifdef _DEBUG
const int theDim[3] = {2, 2, 1};
#else
//const int theDim[3] = {12, 12, 4};
const int theDim[3] = {2, 2, 1};
//const int theDim[3] = {12, 12, 6};
#endif

//const double theCellSize = 0.5;
const double theCellSize = 1.0;

const double fluidDensity = 1.0;
