// Modified by Peter Kutz, 2011 and 2012.

#include "mac_grid.h"
#include "open_gl_headers.h"
#include "camera.h"
#include "custom_output.h"
#include "constants.h"
#include <math.h>
#include <map>
#include <stdio.h>
#undef max
#undef min
#include <fstream>

#include "dprint.h"

// Globals:
MACGrid target;


// NOTE: x -> cols, z -> rows, y -> stacks
MACGrid::RenderMode MACGrid::theRenderMode = SHEETS;
bool MACGrid::theDisplayVel = false;

#define FOR_EACH_CELL \
  for(int k = 0; k < theDim[MACGrid::Z]; k++)  \
    for(int j = 0; j < theDim[MACGrid::Y]; j++) \
      for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
  for(int k = theDim[MACGrid::Z] - 1; k >= 0; k--)  \
    for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
      for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
  for(int k = 0; k < theDim[MACGrid::Z]+1; k++) \
    for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
      for(int i = 0; i < theDim[MACGrid::X]+1; i++) 

#define FOR_EACH_FACE_X \
  for (int k = 0; k < theDim[MACGrid::Z]; k++) \
    for (int j = 0; j < theDim[MACGrid::Y]; j++) \
      for (int i = 0; i < theDim[MACGrid::X]+1; i++)

#define FOR_EACH_FACE_Y \
  for (int k = 0; k < theDim[MACGrid::Z]; k++) \
    for (int j = 0; j < theDim[MACGrid::Y]+1; j++) \
      for (int i = 0; i < theDim[MACGrid::X]; i++)

#define FOR_EACH_FACE_Z \
  for (int k = 0; k < theDim[MACGrid::Z]+1; k++) \
    for (int j = 0; j < theDim[MACGrid::Y]; j++) \
      for (int i = 0; i < theDim[MACGrid::X]; i++)


MACGrid::MACGrid() {
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig) {
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig) {
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   

   return *this;
}

MACGrid::~MACGrid() {
}

void MACGrid::reset() {
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);

   setUpAMatrix();
}

void MACGrid::initialize() {
  reset();
}

void MACGrid::updateSources() {
  static int count = 0;
  // TODO: Set initial values for density, temperature, and velocity.
  // 12x12x1
  /*
  mV(6,1,0) = 1.0;
  mU(7,1,0) = 1.0;
  mU(6,1,0) = -1.0;
  if (count < 100) {
    mD(6,0,0) = 1.0;
  }
  */
  //mV(6,1,0) = 1.0;
  //mU(1,0,0) = 1.0;
  //mT(6,7,0) = 50.0;
  //mD(0,4,0) = 2.0;
  //mU(1,4,0) = 2.0;
  //mV(0,1,0) = 1.0;
  //mT(0,4,0) = 50.0;
  //mD(19,0,0) = 1.0;
  //mU(19,0,0) = -1.0;

  //mT(6,11,0) = -20.0;
  //mD(6,11,0) = 2.0;

  // 2x2x1
  //mV(0,1,0) = 1.0;
  //mU(1,0,0) = 1.0;
  //mT(0,1,0) = 50.0;
  //mD(0,0,0) = 1.0;

  mU(1,0,0) = 2.0;
  //mU(38,39,0) = -2.0;
  mD(0,0,0) = 1.0;
  mD(39,39,0) = 1.0;
  mT(0,0,0) = 50.0;
  mT(39,39,0) = -50.0;

  /*
  mD(4,4,4) = 1.0;
  mD(4,5,4) = 1.0;
  mD(4,4,5) = 1.0;
  mD(4,5,5) = 1.0;
  mT(4,4,4) = 100.0;
  mT(4,5,4) = 100.0;
  mT(4,4,5) = 100.0;
  mT(4,5,5) = 100.0;

  mU(5,4,4) = 4.0;
  mU(5,5,4) = 4.0;
  mU(5,4,5) = 4.0;
  mU(5,5,5) = 4.0;
  */

  if (count >= 1) {
    //exit(0);
  }

  count += 1;
}

void MACGrid::advectVelocity(double dt) {
  // TODO: Calculate new velocities and store in target.
  target.mU = mU;
  target.mV = mV;
  target.mW = mW;

  GridData iX(mU);
  GridData iY(mU);
  GridData iZ(mU);

  // iterate over each Z face
  FOR_EACH_FACE_X {
    // actual grid point
    vec3 pt(i, j, k);
    pt *= theCellSize;
    pt[1] += 0.5 * theCellSize;
    pt[2] += 0.5 * theCellSize;

    // get full dimensional velocity
    vec3 vel = getVelocity(pt);
    //pt.Print("pt: ");
    //vel.Print("vel: ");
    //fflush(stdout);
    iX(i,j,k) = vel[0];
    iY(i,j,k) = vel[1];
    iZ(i,j,k) = vel[2];

    // do backwards euler step
    vec3 bpt = pt - dt * vel;

    // interpolate new velocity
    vec3 nvel = getVelocity(bpt);

    // store new velocity
    target.mU(i, j, k) = nvel[0];
  }


  //printf("y:\n");
  // iterate over each Y face
  FOR_EACH_FACE_Y {
    // actual grid point
    vec3 pt(i, j, k);
    pt *= theCellSize;
    pt[0] += 0.5 * theCellSize;
    pt[2] += 0.5 * theCellSize;

    // get full dimensional velocity
    vec3 vel = getVelocity(pt);
    //pt.Print("pt: ");
    //vel.Print("vel: ");
    //fflush(stdout);

    // do backwards euler step
    vec3 bpt = pt - dt * vel;

    // interpolate new velocity
    vec3 nvel = getVelocity(bpt);

    // store new velocity
    target.mV(i, j, k) = nvel[1];
  }

  // iterate over each Z face
  FOR_EACH_FACE_Z {
    // actual grid point
    vec3 pt(i, j, k);
    pt *= theCellSize;
    pt[0] += 0.5 * theCellSize;
    pt[1] += 0.5 * theCellSize;

    // get full dimensional velocity
    vec3 vel = getVelocity(pt);
    
    // do backwards euler step
    vec3 bpt = pt - dt * vel;

    // interpolate new velocity
    vec3 nvel = getVelocity(bpt);

    // store new velocity
    target.mW(i, j, k) = nvel[2];
  }

  #ifdef __DPRINT__
  #ifdef __DPRINT_ADVVEL__
  printf("***********************************************************************************\n");
  printf("Advected Velocities:\n");
  printf("mU:\n");
  print_grid_data(mU);
  printf("mV:\n");
  print_grid_data(mV);
  printf("mW:\n");
  print_grid_data(mW);
  printf("iX:\n");
  print_grid_data(iX);
  printf("iY:\n");
  print_grid_data(iY);
  printf("iZ:\n");
  print_grid_data(iZ);
  printf("target.mU:\n");
  print_grid_data(target.mU);
  printf("target.mV:\n");
  print_grid_data(target.mV);
  printf("target.mW:\n");
  print_grid_data(target.mW);
  #endif
  #endif

  // Then save the result to our object.
  mU = target.mU;
  mV = target.mV;
  mW = target.mW;
}

void MACGrid::advectTemperature(double dt) {
  // TODO: Calculate new temp and store in target.
  target.mT = mT;

  // temperature is stored per cell
  FOR_EACH_CELL {
    // compute world point
    vec3 pt(i, j, k);
    pt *= theCellSize;
    pt[0] += 0.5 * theCellSize;
    pt[1] += 0.5 * theCellSize;
    pt[2] += 0.5 * theCellSize;

    // get interpolated velocity at the world point
    vec3 vel = getVelocity(pt);

    // euler step to get previous position
    vec3 bpt = pt - dt * vel;

    // get the interpolated temperature
    double nT = getTemperature(bpt);

    // store new temperature
    target.mT(i,j,k) = nT;
  }

  #ifdef __DPRINT__
  #ifdef __DPRINT_ADVTEMP__
  printf("***********************************************************************************\n");
  printf("Advected Temperature:\n");
  printf("mT:\n");
  print_grid_data(mT);
  printf("target.mT:\n");
  print_grid_data(target.mT);
  #endif
  #endif

  // Then save the result to our object.
  mT = target.mT;
}

void MACGrid::advectDensity(double dt) {
  // TODO: Calculate new densitities and store in target.
  target.mD = mD;

  // density is stored per cell
  FOR_EACH_CELL {
    // compute world point
    vec3 pt(i, j, k);
    pt *= theCellSize;
    pt[0] += 0.5 * theCellSize;
    pt[1] += 0.5 * theCellSize;
    pt[2] += 0.5 * theCellSize;

    // get interpolated velocity at the world point
    vec3 vel = getVelocity(pt);

    // euler step to get previous position
    vec3 bpt = pt - dt * vel;

    // get the interpolated density
    double nD = getDensity(bpt);

    // store new density 
    target.mD(i,j,k) = nD;
  }

  #ifdef __DPRINT__
  #ifdef __DPRINT_ADVDENS__
  printf("***********************************************************************************\n");
  printf("Advected Density:\n");
  printf("mD:\n");
  print_grid_data(mD);
  printf("target.mD:\n");
  print_grid_data(target.mD);
  #endif
  #endif

  // Then save the result to our object.
  mD = target.mD;
}

void MACGrid::computeBouyancy(double dt) {
  // TODO: Calculate bouyancy and store in target.
  target.mV = mV;
  double a = buoyAlpha;
  double B = buoyBeta;
  // TODO: what is the mass?
  double mass = 1.0;

  FOR_EACH_FACE_Y {
    if (j == 0 || j == theDim[MACGrid::Y]) {
      // do not update boundary
      continue;
    }
    // get world point for face
    vec3 pt(i,j,k);
    pt *= theCellSize;
    pt[0] += 0.5 * theCellSize;
    pt[2] += 0.5 * theCellSize;

    // get interpolated temperature at face
    double T = getTemperature(pt);
    
    // get interpolated density
    double s = getDensity(pt);

    // buoyancy only affects vertical velocity
    double f = -a * s + B * (T - Tamb);

    // update target velocity
    // TODO: convert mass to acceleration?
    target.mV(i,j,k) = mV(i,j,k) + dt * f;
    //target.mV(i,j,k) = mV(i,j,k) + dt * f/(s*theCellSize*theCellSize + 10e-20);
  }

  #ifdef __DPRINT__
  #ifdef __DPRINT_BUOY__
  printf("***********************************************************************************\n");
  printf("Buoancy:\n");
  printf("mV:\n");
  print_grid_data(mV);
  printf("target.mV:\n");
  print_grid_data(target.mV);
  #endif
  #endif

  // Then save the result to our object.
  mV = target.mV;
}

void MACGrid::computeVorticityConfinement(double dt) {
  // TODO: Calculate vorticity confinement forces.
  // Apply the forces to the current velocity and store the result in target.
  target.mU = mU;
  target.mV = mV;
  target.mW = mW;

  // vorticity confinement coefficient
  double e = vorticityEpsilon;

  // get dimensions
  vec3 dim = mT.getDim();
  double xdim = dim[0];
  double ydim = dim[1];
  double zdim = dim[2];

  GridData wX(mD);
  GridData wY(mD);
  GridData wZ(mD);
  GridData cU(mD);
  GridData cV(mD);
  GridData cW(mD);
  // compute central differences
  FOR_EACH_CELL {
    // get world point of the cell
    vec3 pt(i,j,k);
    pt *= theCellSize;
    pt += vec3(1.0,1.0,1.0)*(0.5*theCellSize);

    // interpolated velocity
    vec3 vel = getVelocity(pt);
    
    // get neighbor positions and velocities
    vec3 ptIplus1 = pt;
    ptIplus1[0] += 1.0 * theCellSize;
    vec3 ptJplus1 = pt;
    ptJplus1[1] += 1.0 * theCellSize;
    vec3 ptKplus1 = pt;
    ptKplus1[2] += 1.0 * theCellSize;
    vec3 ptIminus1 = pt;
    ptIminus1[0] -= 1.0 * theCellSize;
    vec3 ptJminus1 = pt;
    ptJminus1[1] -= 1.0 * theCellSize;
    vec3 ptKminus1 = pt;
    ptKminus1[2] -= 1.0 * theCellSize;


    vec3 velIplus1 = getVelocity(ptIplus1);
    vec3 velJplus1 = getVelocity(ptJplus1);
    vec3 velKplus1 = getVelocity(ptKplus1);
    vec3 velIminus1 = getVelocity(ptIminus1);
    vec3 velJminus1 = getVelocity(ptJminus1);
    vec3 velKminus1 = getVelocity(ptKminus1);

    // compute w (voricity)
    // w_{i,j,k} = (1/(2*dx)) * [ (w_{i,j+1,k} - w_{i,j-1,k}) - (v_{i,j,k+1} - v_{i,j,k-1}),
    //                            (u_{i,j,k+1} - u_{i,j,k-1}) - (w_{i+1,j,k} - w_{i-1,j,k}),
    //                            (v_{i+1,j,k} - v_{i-1,j,k}) - (u_{i,j+1,k} - u_{i,j-1,k}) ];
    wX(i,j,k) = ((velJplus1[2] - velJminus1[2]) - (velKplus1[1] - velKminus1[1])) / (2 * theCellSize);
    wY(i,j,k) = ((velKplus1[0] - velKminus1[0]) - (velIplus1[2] - velIminus1[2])) / (2 * theCellSize);
    wZ(i,j,k) = ((velIplus1[1] - velIminus1[1]) - (velJplus1[0] - velJminus1[0])) / (2 * theCellSize);
  }

  #ifdef __DPRINT__
  #ifdef __DPRINT_VORT__
  printf("***********************************************************************************\n");
  printf("Vortices:\n");
  printf("mU:\n");
  print_grid_data(mU);
  printf("mV:\n");
  print_grid_data(mV);
  printf("mW:\n");
  print_grid_data(mW);
  printf("wX:\n");
  print_grid_data(wX);
  printf("wY:\n");
  print_grid_data(wY);
  printf("wZ:\n");
  print_grid_data(wZ);
  #endif
  #endif

  GridData gwX(mD);
  GridData gwY(mD);
  GridData gwZ(mD);
  GridData fX(mD);
  GridData fY(mD);
  GridData fZ(mD);
  GridDataX dmU(mU);
  GridDataY dmV(mV);
  GridDataZ dmW(mW);
  // compute gradient of w
  FOR_EACH_CELL {
    // TODO: what to do about difference that are outside the grid?
    //  gradW_{i,j,k} = ( (|w_{i+1,j,k}| - |w_{i-1,j,k}|), (|w_{i,j+1,k}| - |w_{i,j-1,k}|), (|w_{i,j,k+1}| - |w_{i,j,k-1}|) )/(2*dx)
    //    where |x| = 1-norm
    vec3 w(wX(i,j,k), wY(i,j,k), wZ(i,j,k));
    
    // get neighboring cells
    //  if border use, replicate the cell data
    vec3 wIplus1(w), wJplus1(w), wKplus1(w);
    vec3 wIminus1(w), wJminus1(w), wKminus1(w);
    if (i + 1 < xdim) {
      wIplus1 = vec3(wX(i+1,j,k), wY(i+1,j,k), wZ(i+1,j,k));
    } 
    if (j + 1 < ydim) {
      wJplus1 = vec3(wX(i,j+1,k), wY(i,j+1,k), wZ(i,j+1,k));
    }
    if (k + 1 < zdim) {
      wKplus1 = vec3(wX(i,j,k+1), wY(i,j,k+1), wZ(i,j,k+1));
    }
    if (i - 1 >= 0) {
      wIminus1 = vec3(wX(i-1,j,k), wY(i-1,j,k), wZ(i-1,j,k));
    }
    if (j - 1 >= 0) {
      wJminus1 = vec3(wX(i,j-1,k), wY(i,j-1,k), wZ(i,j-1,k));
    }
    if (k - 1 >= 0) {
      wKminus1 = vec3(wX(i,j,k-1), wY(i,j,k-1), wZ(i,j,k-1));
    }

    /*
    gwX(i,j,k) = (wIplus1.SqrLength() - wIminus1.SqrLength())/(2*theCellSize);
    gwY(i,j,k) = (wJplus1.SqrLength() - wJminus1.SqrLength())/(2*theCellSize);
    gwZ(i,j,k) = (wKplus1.SqrLength() - wKminus1.SqrLength())/(2*theCellSize);
    */

    

    gwX(i,j,k) = ((abs(wIplus1[0]) + abs(wIplus1[1]) + abs(wIplus1[2])) - (abs(wIminus1[0]) + abs(wIminus1[1]) + abs(wIminus1[2]))) / (2*theCellSize); 
    gwY(i,j,k) = ((abs(wJplus1[0]) + abs(wJplus1[1]) + abs(wJplus1[2])) - (abs(wJminus1[0]) + abs(wJminus1[1]) + abs(wJminus1[2]))) / (2*theCellSize); 
    gwZ(i,j,k) = ((abs(wKplus1[0]) + abs(wKplus1[1]) + abs(wKplus1[2])) - (abs(wKminus1[0]) + abs(wKminus1[1]) + abs(wKminus1[2]))) / (2*theCellSize); 
    


    
    // normalize
    vec3 gw(gwX(i,j,k), gwY(i,j,k), gwZ(i,j,k));
    double div = gw.Length() + 10e-20;
    gwX(i,j,k) = gwX(i,j,k) / div;
    gwY(i,j,k) = gwY(i,j,k) / div;
    gwZ(i,j,k) = gwZ(i,j,k) / div;

    // compute force
    // f = epsilon * dx * (N cross w);
    vec3 n(gwX(i,j,k), gwY(i,j,k), gwZ(i,j,k));
    vec3 fconf = vorticityEpsilon * (n ^ w);

    // store forces
    fX(i,j,k) = fconf[0];
    fY(i,j,k) = fconf[1];
    fZ(i,j,k) = fconf[2];
  }

  // update velocities
  FOR_EACH_FACE_X {
    // do not add forces to the edge faces
    if (i > 0 && i < xdim) {
      // get world pt
      vec3 pt(i,j,k);
      pt *= theCellSize;
      pt += vec3(0.0,1.0,1.0) * (0.5*theCellSize);
      
      // get interpolated force at face
      double f = fX.interpolate(pt);

      // update velocity
      target.mU(i,j,k) = mU(i,j,k) + dt * f;
      dmU(i,j,k) = dt * f;
    } else {
      dmU(i,j,k) = 0;
    }
  }

  FOR_EACH_FACE_Y {
    // do not add forces to the edge faces
    if (j > 0 && j < ydim) {
      // get world pt
      vec3 pt(i,j,k);
      pt *= theCellSize;
      pt += vec3(1.0,0.0,1.0) * (0.5*theCellSize);
      
      // get interpolated force at face
      double f = fY.interpolate(pt);

      // update velocity
      target.mV(i,j,k) = mV(i,j,k) + dt * f;
      dmV(i,j,k) = dt * f;
    } else {
      dmV(i,j,k) = 0;
    }
  }

  FOR_EACH_FACE_Z {
    // do not add forces to the edge faces
    if (k > 0 && k < zdim) {
      // get world pt
      vec3 pt(i,j,k);
      pt *= theCellSize;
      pt += vec3(1.0,1.0,0.0) * (0.5*theCellSize);
      
      // get interpolated force at face
      double f = fZ.interpolate(pt);

      // update velocity
      target.mW(i,j,k) = mW(i,j,k) + dt * f;
      dmW(i,j,k) = dt * f;
    } else {
      dmW(i,j,k) = 0;
    }
  }

  #ifdef __DPRINT__
  #ifdef __DPRINT_VORT__
  printf("-----------------------------------------------------\n");
  printf("gwX:\n");
  print_grid_data(gwX);
  printf("gwY:\n");
  print_grid_data(gwY);
  printf("gwZ:\n");
  print_grid_data(gwZ);
  printf("fX:\n");
  print_grid_data(fX);
  printf("fY:\n");
  print_grid_data(fY);
  printf("fZ:\n");
  print_grid_data(fZ);
  /*
  printf("dmU:\n");
  print_grid_data(dmU);
  printf("dmV:\n");
  print_grid_data(dmV);
  printf("dmW:\n");
  print_grid_data(dmW);
  */
  #endif
  #endif


  // Then save the result to our object.
  mU = target.mU;
  mV = target.mV;
  mW = target.mW;
}

void MACGrid::addExternalForces(double dt) {
   computeBouyancy(dt);
   computeVorticityConfinement(dt);
}

void MACGrid::project(double dt) {
  // TODO: Solve Ap = d for pressure.
  // 1. Construct d
  // 2. Construct A
  // 3. Solve for p
  // Subtract pressure from our velocity and save in target.
  target.mP = mP;
  target.mU = mU;
  target.mV = mV;
  target.mW = mW;


  // construct d (RHS)
  //  change in velocity
  //  d(i,j,k) = -((rho)*(dx^2)/dt) * (change in total velocities)
  GridData d = mP;
  FOR_EACH_CELL {
    // compute constant
    // TODO: rho is just the denisty right?
    //double c = -1.0 * mD(i,j,k) * theCellSize * theCellSize / dt;
    double c = -1.0 * fluidDensity * theCellSize * theCellSize / dt;

    // compute change in x velocity
    double u = (mU(i+1,j,k) - mU(i,j,k)) / theCellSize;
    // compute change in y velocity
    double v = (mV(i,j+1,k) - mV(i,j,k)) / theCellSize;
    // compute change in y velocity
    double w = (mW(i,j,k+1) - mW(i,j,k)) / theCellSize;

    // store sum in d
    d(i,j,k) = c * (u + v + w);
  }


  // construct A
  //  matrix consisting of pressure neighbor information
  //  already constructed for boundaries
  //  TODO: add support for obstacles in the grid
  GridDataMatrix A = AMatrix;


  // solve for new pressures such that the fluid remains incompressible
  // TODO: what maxIterations and tolerance to use?
  bool ret = conjugateGradient(A, target.mP, d, 100000, 0.00001);


  #ifdef __DPRINT__
  #ifdef __DPRINT_PROJECT__
  printf("***********************************************************************************\n");
  printf("Project: %s\n", ret ? "true" : "false");
  //printf("A:\n");
  //A.print();
  printf("d:\n");
  print_grid_data_as_column(d);
  printf("mP:\n");
  print_grid_data_as_column(target.mP);
  #endif
  #endif


  // update velocities from new pressures
  // vn = v - dt * (1/rho) * dP
  FOR_EACH_CELL {
    // update velocity faces (+1 cell index)
    //  boundaries should not change

    // update x velocity
    if (i < theDim[MACGrid::X] - 1) {
      target.mU(i+1,j,k) = mU(i+1,j,k) - dt * (target.mP(i+1,j,k) - target.mP(i,j,k));
    }

    // update y velocity
    if (j < theDim[MACGrid::Y] - 1) {
      target.mV(i,j+1,k) = mV(i,j+1,k) - dt * (target.mP(i,j+1,k) - target.mP(i,j,k));
    }
    
    // update z velocity
    if (k < theDim[MACGrid::Z] - 1) {
      target.mW(i,j,k+1) = mW(i,j,k+1) - dt * (target.mP(i,j,k+1) - target.mP(i,j,k));
    }
  }


  // Then save the result to our object
  mP = target.mP;
  mU = target.mU;
  mV = target.mV;
  mW = target.mW;
  // IMPLEMENT THIS AS A SANITY CHECK: assert (checkDivergence());
  //assert(checkDivergence());
  if (!checkDivergence()) {
    //exit(0);
  }
}

bool MACGrid::checkDivergence() {
  double sum = 0;
  bool ret = true;
  FOR_EACH_FACE_X {
    sum += mU(i,j,k);
  }
  if (abs(sum) > 10e-6) {
    //printf("X vel DIVERGENT: %f\n", sum); 
    ret = false;
  }
  FOR_EACH_FACE_Y {
    sum += mV(i,j,k);
  }
  if (abs(sum) > 10e-6) {
    //printf("Y vel DIVERGENT: %f\n", sum); 
    ret = false;
  }
  FOR_EACH_FACE_Z {
    sum += mW(i,j,k);
  }
  if (abs(sum) > 10e-6) {
    //printf("Z vel DIVERGENT: %f\n", sum); 
    ret = false;
  }
  fflush(stdout);
  return ret; 
}


vec3 MACGrid::getVelocity(const vec3& pt) {
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt) {
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt) {
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt) {
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt) {
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt) {
   return mD.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k) {
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}

bool MACGrid::isValidCell(int i, int j, int k) {
  if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
    return false;
  }

  if (i < 0 || j < 0 || k < 0) {
    return false;
  }

  return true;
}

void MACGrid::setUpAMatrix() {

  FOR_EACH_CELL {

    int numFluidNeighbors = 0;
    if (i-1 >= 0) {
      AMatrix.plusI(i-1,j,k) = -1;
      numFluidNeighbors++;
    }
    if (i+1 < theDim[MACGrid::X]) {
      AMatrix.plusI(i,j,k) = -1;
      numFluidNeighbors++;
    }
    if (j-1 >= 0) {
      AMatrix.plusJ(i,j-1,k) = -1;
      numFluidNeighbors++;
    }
    if (j+1 < theDim[MACGrid::Y]) {
      AMatrix.plusJ(i,j,k) = -1;
      numFluidNeighbors++;
    }
    if (k-1 >= 0) {
      AMatrix.plusK(i,j,k-1) = -1;
      numFluidNeighbors++;
    }
    if (k+1 < theDim[MACGrid::Z]) {
      AMatrix.plusK(i,j,k) = -1;
      numFluidNeighbors++;
    }
    // Set the diagonal:
    AMatrix.diag(i,j,k) = numFluidNeighbors;
  }
}







/////////////////////////////////////////////////////////////////////

bool MACGrid::conjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
  // Solves Ap = d for p.

  FOR_EACH_CELL {
    p(i,j,k) = 0.0; // Initial guess p = 0. 
  }

  GridData r = d; // Residual vector.

  GridData z; z.initialize();
  // TODO: Apply a preconditioner here.
  // For now, just bypass the preconditioner:
  z = r;

  GridData s = z; // Search vector;

  double sigma = dotProduct(z, r);

  for (int iteration = 0; iteration < maxIterations; iteration++) {

    double rho = sigma;

    apply(A, s, z);

    double alpha = rho/dotProduct(z, s);

    GridData alphaTimesS; alphaTimesS.initialize();
    multiply(alpha, s, alphaTimesS);
    add(p, alphaTimesS, p);

    GridData alphaTimesZ; alphaTimesZ.initialize();
    multiply(alpha, z, alphaTimesZ);
    subtract(r, alphaTimesZ, r);

    if (maxMagnitude(r) <= tolerance) {
      //PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
      return true;
    }

    // TODO: Apply a preconditioner here.
    // For now, just bypass the preconditioner:
    z = r;

    double sigmaNew = dotProduct(z, r);

    double beta = sigmaNew / rho;

    GridData betaTimesS; betaTimesS.initialize();
    multiply(beta, s, betaTimesS);
    add(z, betaTimesS, s);
    //s = z + beta * s;

    sigma = sigmaNew;
  }

  PRINT_LINE( "PCG didn't converge!" );
  return false;

}

double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
  
  double result = 0.0;

  FOR_EACH_CELL {
    result += vector1(i,j,k) * vector2(i,j,k);
  }

  return result;
}

void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
  
  FOR_EACH_CELL {
    result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
  }

}

void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
  
  FOR_EACH_CELL {
    result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
  }

}

void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
  
  FOR_EACH_CELL {
    result(i,j,k) = scalar * vector(i,j,k);
  }

}

double MACGrid::maxMagnitude(const GridData & vector) {
  
  double result = 0.0;

  FOR_EACH_CELL {
    if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
  }

  return result;
}

void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
  
  FOR_EACH_CELL { // For each row of the matrix.

    double diag = 0;
    double plusI = 0;
    double plusJ = 0;
    double plusK = 0;
    double minusI = 0;
    double minusJ = 0;
    double minusK = 0;

    diag = matrix.diag(i,j,k) * vector(i,j,k);
    if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
    if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
    if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
    if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
    if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
    if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

    result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
  }

}




/////////////////////////////////////////////////////////////////////

void MACGrid::saveSmoke(const char* fileName) {
  std::ofstream fileOut(fileName);
  if (fileOut.is_open()) {
    FOR_EACH_CELL {
      fileOut << mD(i,j,k) << std::endl;
    }
    fileOut.close();
  }
}





/////////////////////////////////////////////////////////////////////

void MACGrid::draw(const Camera& c) {   
   drawWireGrid();
   if (theDisplayVel) drawVelocities();   
   if (theRenderMode == CUBES) drawSmokeCubes(c);
   else drawSmoke(c);
}

void MACGrid::drawVelocities() {
  // Draw line at each center
  //glColor4f(0.0, 1.0, 0.0, 1.0); // Use this if you want the lines to be a single color.
  glBegin(GL_LINES);
    FOR_EACH_CELL {
      vec3 pos = getCenter(i,j,k);
      vec3 vel = getVelocity(pos);
      if (vel.Length() > 0.0001) {
        // Un-comment the line below if you want all of the velocity lines to be the same length.
        //vel.Normalize();
        vel *= theCellSize/2.0;
        vel += pos;
        glColor4f(1.0, 1.0, 0.0, 1.0);
        glVertex3dv(pos.n);
        glColor4f(0.0, 1.0, 0.0, 1.0);
        glVertex3dv(vel.n);
      }
    }
  glEnd();
}

vec4 MACGrid::getRenderColor(int i, int j, int k) {
  // Modify this if you want to change the smoke color, or modify it based on other smoke properties.
  double value = mD(i, j, k); 
  return vec4(1.0, 1.0, 1.0, value);
}

vec4 MACGrid::getRenderColor(const vec3& pt) {
  // TODO: Modify this if you want to change the smoke color, or modify it based on other smoke properties.
  double value = getDensity(pt); 

  return vec4(1.0, 1.0, 1.0, value);

  /*
  double temp = getTemperature(pt);
  temp = max(Tmin, min(temp, Tmax));
  double v = 0.75;
  double s = 0.75;

  double hh = 6 * (temp - Tmin)/(Tmax - Tmin);
  int i = hh;
  double ff = hh - i;
  double p = v * (1.0 - s);
  double q = v * (1.0 - (s * ff));
  double t = v * (1.0 - (s * (1.0 - ff)));

  double r = 0;
  double g = 0;
  double b = 0;

  switch (i) { 
    case 0:
      r = v;
      g = t;
      b = p;
    case 1:
      r = q;
      g = v;
      b = p;
    case 2:
      r = p;
      g = v;
      b = t;
    case 3:
      r = p;
      g = q;
      b = v;
    case 4:
      r = t;
      g = p;
      b = v;
    case 5:
      r = v;
      g = p;
      b = q;
  }

  //return vec4(r, g, b, value);
  return vec4(b, g, r, value);
  */
}

void MACGrid::drawZSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double startk = back - stepsize;
   double endk = 0;
   double stepk = -theCellSize;

   if (!backToFront)
   {
      startk = 0;
      endk = back;   
      stepk = theCellSize;
   }

   for (double k = startk; backToFront? k > endk : k < endk; k += stepk)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double i = 0.0; i <= right; i += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double i = right; i >= 0.0; i -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawXSheets(bool backToFront)
{
   // Draw K Sheets from back to front
   double back =  (theDim[2])*theCellSize;
   double top  =  (theDim[1])*theCellSize;
   double right = (theDim[0])*theCellSize;
  
   double stepsize = theCellSize*0.25;

   double starti = right - stepsize;
   double endi = 0;
   double stepi = -theCellSize;

   if (!backToFront)
   {
      starti = 0;
      endi = right;   
      stepi = theCellSize;
   }

   for (double i = starti; backToFront? i > endi : i < endi; i += stepi)
   {
     for (double j = 0.0; j < top; )
      {
         glBegin(GL_QUAD_STRIP);
         for (double k = 0.0; k <= back; k += stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;

         glBegin(GL_QUAD_STRIP);
         for (double k = back; k >= 0.0; k -= stepsize)
         {
            vec3 pos1 = vec3(i,j,k); 
            vec3 pos2 = vec3(i, j+stepsize, k); 

            vec4 color1 = getRenderColor(pos1);
            vec4 color2 = getRenderColor(pos2);

            glColor4dv(color1.n);
            glVertex3dv(pos1.n);

            glColor4dv(color2.n);
            glVertex3dv(pos2.n);
         } 
         glEnd();
         j+=stepsize;
      }
   }
}

void MACGrid::drawSmoke(const Camera& c)
{
   vec3 eyeDir = c.getBackward();
   double zresult = fabs(Dot(eyeDir, vec3(1,0,0)));
   double xresult = fabs(Dot(eyeDir, vec3(0,0,1)));
   //double yresult = fabs(Dot(eyeDir, vec3(0,1,0)));

   if (zresult < xresult)
   {      
      drawZSheets(c.getPosition()[2] < 0);
   }
   else 
   {
      drawXSheets(c.getPosition()[0] < 0);
   }
}

void MACGrid::drawSmokeCubes(const Camera& c)
{
   std::multimap<double, MACGrid::Cube, std::greater<double> > cubes;
   FOR_EACH_CELL
   {
      MACGrid::Cube cube;
      cube.color = getRenderColor(i,j,k);
      cube.pos = getCenter(i,j,k);
      cube.dist = DistanceSqr(cube.pos, c.getPosition());
      cubes.insert(make_pair(cube.dist, cube));
   } 

   // Draw cubes from back to front
   std::multimap<double, MACGrid::Cube, std::greater<double> >::const_iterator it;
   for (it = cubes.begin(); it != cubes.end(); ++it)
   {
      drawCube(it->second);
   }
}

void MACGrid::drawWireGrid()
{
   // Display grid in light grey, draw top & bottom

   double xstart = 0.0;
   double ystart = 0.0;
   double zstart = 0.0;
   double xend = theDim[0]*theCellSize;
   double yend = theDim[1]*theCellSize;
   double zend = theDim[2]*theCellSize;

   glPushAttrib(GL_LIGHTING_BIT | GL_LINE_BIT);
      glDisable(GL_LIGHTING);
      glColor3f(0.25, 0.25, 0.25);

      glBegin(GL_LINES);
      for (int i = 0; i <= theDim[0]; i++)
      {
         double x = xstart + i*theCellSize;
         glVertex3d(x, ystart, zstart);
         glVertex3d(x, ystart, zend);

         glVertex3d(x, yend, zstart);
         glVertex3d(x, yend, zend);
      }

      for (int i = 0; i <= theDim[2]; i++)
      {
         double z = zstart + i*theCellSize;
         glVertex3d(xstart, ystart, z);
         glVertex3d(xend, ystart, z);

         glVertex3d(xstart, yend, z);
         glVertex3d(xend, yend, z);
      }

      glVertex3d(xstart, ystart, zstart);
      glVertex3d(xstart, yend, zstart);

      glVertex3d(xend, ystart, zstart);
      glVertex3d(xend, yend, zstart);

      glVertex3d(xstart, ystart, zend);
      glVertex3d(xstart, yend, zend);

      glVertex3d(xend, ystart, zend);
      glVertex3d(xend, yend, zend);
      glEnd();
   glPopAttrib();

   glEnd();
}

#define LEN 0.5
void MACGrid::drawFace(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);
      glEnd();
   glPopMatrix();
}

void MACGrid::drawCube(const MACGrid::Cube& cube)
{
   glColor4dv(cube.color.n);
   glPushMatrix();
      glTranslated(cube.pos[0], cube.pos[1], cube.pos[2]);      
      glScaled(theCellSize, theCellSize, theCellSize);
      glBegin(GL_QUADS);
         glNormal3d( 0.0, -1.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN,  LEN);
         glVertex3d( LEN, -LEN, -LEN);         

         glNormal3d( 0.0,  0.0, -0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN,  LEN, -LEN);
         glVertex3d( LEN,  LEN, -LEN);
         glVertex3d( LEN, -LEN, -LEN);

         glNormal3d(-1.0,  0.0,  0.0);
         glVertex3d(-LEN, -LEN, -LEN);
         glVertex3d(-LEN, -LEN,  LEN);
         glVertex3d(-LEN,  LEN,  LEN);
         glVertex3d(-LEN,  LEN, -LEN);

         glNormal3d( 0.0, 1.0,  0.0);
         glVertex3d(-LEN, LEN, -LEN);
         glVertex3d(-LEN, LEN,  LEN);
         glVertex3d( LEN, LEN,  LEN);
         glVertex3d( LEN, LEN, -LEN);

         glNormal3d( 0.0,  0.0, 1.0);
         glVertex3d(-LEN, -LEN, LEN);
         glVertex3d(-LEN,  LEN, LEN);
         glVertex3d( LEN,  LEN, LEN);
         glVertex3d( LEN, -LEN, LEN);

         glNormal3d(1.0,  0.0,  0.0);
         glVertex3d(LEN, -LEN, -LEN);
         glVertex3d(LEN, -LEN,  LEN);
         glVertex3d(LEN,  LEN,  LEN);
         glVertex3d(LEN,  LEN, -LEN);
      glEnd();
   glPopMatrix();
}
