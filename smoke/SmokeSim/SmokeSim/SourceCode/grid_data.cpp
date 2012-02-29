#include "grid_data.h"

//#define __CUBIC_INTERP__
#define __BRIDSON_NOTATION__

GridData::GridData() :
   mDfltValue(0.0), mMax(0.0,0.0,0.0)
{
}

GridData::GridData(const GridData& orig) :
   mDfltValue(orig.mDfltValue)
{
   mData = orig.mData;
   mMax = orig.mMax;
}

GridData::~GridData() 
{
}

std::vector<double>& GridData::data()
{
   return mData;
}

vec3 GridData::getDim()
{
   return vec3(mMax);
}

GridData& GridData::operator=(const GridData& orig)
{
   if (this == &orig)
   {
      return *this;
   }
   mDfltValue = orig.mDfltValue;
   mData = orig.mData;
   mMax = orig.mMax;
   return *this;
}

void GridData::initialize(double dfltValue)
{
   mDfltValue = dfltValue;
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridData::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}

const double GridData::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // HACK: Protect against setting the default value

   if (i< 0 || j<0 || k<0 || 
       i > theDim[0]-1 || 
       j > theDim[1]-1 || 
       k > theDim[2]-1) return dflt;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];

   return mData[col+row+stack];
}

void GridData::getCell(const vec3& pt, int& i, int& j, int& k)
{
   vec3 pos = worldToSelf(pt); 
   i = (int) (pos[0]/theCellSize);
   j = (int) (pos[1]/theCellSize);
   k = (int) (pos[2]/theCellSize);   
}

double GridData::cubic_interp(double fm1, double f0, double f1, double f2, double dfrac) {
  #ifndef __BRIDSON_NOTATION__
  // t - position; f - value
  // f = a_3 * (t - t_k)^3 + a_2 * (t - t_k)^2 + a_1 * (t - tk) + a_0;
  // a_0 = f_k
  // a_1 = d_k
  // a_2 = 3*deltak - 2*d_k - d_{k+1}
  // a_3 = d_k + d_{k+1} - 2*deltak
  // d_k = (f_{k+1} - f_{k-1}) / 2.0
  // deltak = f_{k+1} - f_k


  double deltak = f1 - f0;
  double d_k = (f1 - fm1)/2.0;
  double d_k1 = (f2 - f0)/2.0;
  
  // if sign of d_k or d_{k+1} are different from deltak set them to 0, repectively
  if (SIGN(d_k) != SIGN(deltak)) {
    d_k = 0.0;
  }
  if (SIGN(d_k1) != SIGN(deltak)) {
    d_k1 = 0.0;
  }

  double a_3 = d_k + d_k1 - 2*deltak;
  double a_2 = 3 * deltak - 2 * d_k - d_k1;
  double a_1 = d_k;
  double a_0 = f0;
  
  double f = a_3 * pow(dfrac, 3) + a_2 * pow(dfrac, 2) + a_1 * dfrac + a_0;

  return f;

  #else

  // bridson notation
  // f = f0 + d_0 * x + (3 * deltaf - 2 * d_0 - d_1) * x^2 + (-2 * deltaf + d_0 + d_1) * x^3
  // d_i = (f_{i+1} - f_{i-1})/2.0
  // deltaf = f_{i+1} - f_i

  double d_0 = (f1 - fm1)/2.0;
  double d_1 = (f2 - f0)/2.0;
  double deltaf = f1 - f0;

  if (SIGN(d_0) != SIGN(deltaf)) {
    d_0 = 0.0;
  }
  if (SIGN(d_1) != SIGN(deltaf)) {
    d_1 = 0.0;
  }

  double f = f0 + d_0 * dfrac + (3*deltaf - 2*d_0 - d_1) * pow(dfrac,2) + (-2*deltaf + d_0 + d_1) * pow(dfrac,3);

  return f;
  #endif
}

double GridData::interpolate(const vec3& pt) {

	// TODO: Implement sharper cubic interpolation here.

	// LINEAR INTERPOLATION:
	vec3 pos = worldToSelf(pt);

	int i = (int) (pos[0]/theCellSize);
	int j = (int) (pos[1]/theCellSize);
	int k = (int) (pos[2]/theCellSize);

	double scale = 1.0/theCellSize;  
	double fractx = scale*(pos[0] - i*theCellSize);
	double fracty = scale*(pos[1] - j*theCellSize);
	double fractz = scale*(pos[2] - k*theCellSize);

	assert (fractx < 1.0 && fractx >= 0);
	assert (fracty < 1.0 && fracty >= 0);
	assert (fractz < 1.0 && fractz >= 0);

  #ifdef __CUBIC_INTERP__

  // 8 conn neighborhood

  // want interpolated value of t which is in the interval [t_k:t_{k+1})
  // t - position; f - value
  // f = a_3 * (t - t_k)^3 + a_2 * (t - t_k)^2 + a_1 * (t - tk) + a_0;
  // a_0 = f_k
  // a_1 = d_k
  // a_2 = 3*deltak - 2*d_k - d_{k+1}
  // a_3 = d_k + d_{k+1} - deltak
  // d_k = (f_{k+1} - f_{k-1}) / 2.0
  // deltak = f_{k+1} - f_k


  // box i x,y,z used in simulation display (x-left, y-up, z-in)
  //  centered at (i,j,k)
  double xfrac = fractx;
  double yfrac = fracty;
  double zfrac = fractz;

  ////// Y interpolation
  //// top level
  // center-right column
  double fm1m10 = (*this)(i-1,  j-1,  k);
  double fm100  = (*this)(i-1,    j,  k);
  double fm110  = (*this)(i-1,  j+1,  k);
  double fm120  = (*this)(i-1,  j+2,  k);
  double yim100 = cubic_interp(fm1m10, fm100, fm110, fm120, yfrac);
  // center column
  double f0m10  = (*this)(i,  j-1,  k);
  double f000   = (*this)(i,    j,  k);
  double f010   = (*this)(i,  j+1,  k);
  double f020   = (*this)(i,  j+2,  k);
  double yi000  = cubic_interp(f0m10, f000, f010, f020, yfrac);
  // center-left column
  double f1m10  = (*this)(i+1,  j-1,  k);
  double f100   = (*this)(i+1,    j,  k);
  double f110   = (*this)(i+1,  j+1,  k);
  double f120   = (*this)(i+1,  j+2,  k);
  double yi100  = cubic_interp(f1m10, f100, f110, f120, yfrac);
  // center-2left column
  double f2m10  = (*this)(i+2,  j-1,  k);
  double f200   = (*this)(i+2,    j,  k);
  double f210   = (*this)(i+2,  j+1,  k);
  double f220   = (*this)(i+2,  j+2,  k);
  double yi200  = cubic_interp(f2m10, f200, f210, f220, yfrac);

  // far-right column
  double fm1m11 = (*this)(i-1,  j-1,  k+1);
  double fm101  = (*this)(i-1,    j,  k+1);
  double fm111  = (*this)(i-1,  j+1,  k+1);
  double fm121  = (*this)(i-1,  j+2,  k+1);
  double yim101 = cubic_interp(fm1m11, fm101, fm111, fm121, yfrac);
  // far-center column
  double f0m11  = (*this)(i,  j-1,  k+1);
  double f001   = (*this)(i,    j,  k+1);
  double f011   = (*this)(i,  j+1,  k+1);
  double f021   = (*this)(i,  j+2,  k+1);
  double yi001  = cubic_interp(f0m11, f001, f011, f021, yfrac);
  // far-left column
  double f1m11  = (*this)(i+1,  j-1,  k+1);
  double f101   = (*this)(i+1,    j,  k+1);
  double f111   = (*this)(i+1,  j+1,  k+1);
  double f121   = (*this)(i+1,  j+2,  k+1);
  double yi101  = cubic_interp(f1m11, f101, f111, f121, yfrac);
  // far-2left column
  double f2m11  = (*this)(i+2,  j-1,  k+1);
  double f201   = (*this)(i+2,    j,  k+1);
  double f211   = (*this)(i+2,  j+1,  k+1);
  double f221   = (*this)(i+2,  j+2,  k+1);
  double yi201  = cubic_interp(f2m11, f201, f211, f221, yfrac);

  // near-right column
  double fm1m1m1  = (*this)(i-1,  j-1,  k-1);
  double fm10m1   = (*this)(i-1,    j,  k-1);
  double fm11m1   = (*this)(i-1,  j+1,  k-1);
  double fm12m1   = (*this)(i-1,  j+2,  k-1);
  double yim10m1  = cubic_interp(fm1m1m1, fm10m1, fm11m1, fm12m1, yfrac);
  // near-center column
  double f0m1m1 = (*this)(i,  j-1,  k-1);
  double f00m1  = (*this)(i,    j,  k-1);
  double f01m1  = (*this)(i,  j+1,  k-1);
  double f02m1  = (*this)(i,  j+2,  k-1);
  double yi00m1 = cubic_interp(f0m1m1, f00m1, f01m1, f02m1, yfrac);
  // near-left column
  double f1m1m1 = (*this)(i+1,  j-1,  k-1);
  double f10m1  = (*this)(i+1,    j,  k-1);
  double f11m1  = (*this)(i+1,  j+1,  k-1);
  double f12m1  = (*this)(i+1,  j+2,  k-1);
  double yi10m1 = cubic_interp(f1m1m1, f10m1, f11m1, f12m1, yfrac);
  // near-2left column
  double f2m1m1 = (*this)(i+2,  j-1,  k-1);
  double f20m1  = (*this)(i+2,    j,  k-1);
  double f21m1  = (*this)(i+2,  j+1,  k-1);
  double f22m1  = (*this)(i+2,  j+2,  k-1);
  double yi20m1 = cubic_interp(f2m1m1, f20m1, f21m1, f22m1, yfrac);

  // +2 far-right column
  double fm1m12   = (*this)(i-1,  j-1,  k+2);
  double fm102    = (*this)(i-1,    j,  k+2);
  double fm112    = (*this)(i-1,  j+1,  k+2);
  double fm122    = (*this)(i-1,  j+2,  k+2);
  double yim102   = cubic_interp(fm1m12, fm102, fm112, fm122, yfrac);
  // +2 far-center column
  double f0m12    = (*this)(i,  j-1,  k+2);
  double f002     = (*this)(i,    j,  k+2);
  double f012     = (*this)(i,  j+1,  k+2);
  double f022     = (*this)(i,  j+2,  k+2);
  double yi002    = cubic_interp(f0m12, f002, f012, f022, yfrac);
  // +2 far-left column
  double f1m12    = (*this)(i+1,  j-1,  k+2);
  double f102     = (*this)(i+1,    j,  k+2);
  double f112     = (*this)(i+1,  j+1,  k+2);
  double f122     = (*this)(i+1,  j+2,  k+2);
  double yi102    = cubic_interp(f1m12, f102, f112, f122, yfrac);
  // +2 far-2left column
  double f2m12    = (*this)(i+2,  j-1,  k+2);
  double f202     = (*this)(i+2,    j,  k+2);
  double f212     = (*this)(i+2,  j+1,  k+2);
  double f222     = (*this)(i+2,  j+2,  k+2);
  double yi202    = cubic_interp(f2m12, f202, f212, f222, yfrac);

  //// Z interpolation
  // right column
  double zim100 = cubic_interp(yim10m1, yim100, yim101, yim102, zfrac);
  // center column
  double zi000 = cubic_interp(yi00m1, yi000, yi001, yi002, zfrac);
  // left column
  double zi100 = cubic_interp(yi10m1, yi100, yi101, yi102, zfrac);
  // +2 left column
  double zi200 = cubic_interp(yi20m1, yi200, yi201, yi202, zfrac);

  //// X interpolation
  double xi000 = cubic_interp(zim100, zi000, zi100, zi200, xfrac);

  return xi000;
    
  #else
	// Y @ low X, low Z:
	double tmp1 = (*this)(i,j,k);
	double tmp2 = (*this)(i,j+1,k);
	// Y @ high X, low Z:
	double tmp3 = (*this)(i+1,j,k);
	double tmp4 = (*this)(i+1,j+1,k);

	// Y @ low X, high Z:
	double tmp5 = (*this)(i,j,k+1);
	double tmp6 = (*this)(i,j+1,k+1);
	// Y @ high X, high Z:
	double tmp7 = (*this)(i+1,j,k+1);
	double tmp8 = (*this)(i+1,j+1,k+1);

	// Y @ low X, low Z
	double tmp12 = LERP(tmp1, tmp2, fracty);
	// Y @ high X, low Z
	double tmp34 = LERP(tmp3, tmp4, fracty);

	// Y @ low X, high Z
	double tmp56 = LERP(tmp5, tmp6, fracty);
	// Y @ high X, high Z
	double tmp78 = LERP(tmp7, tmp8, fracty);

	// X @ low Z
	double tmp1234 = LERP (tmp12, tmp34, fractx);
	// X @ high Z
	double tmp5678 = LERP (tmp56, tmp78, fractx);

	// Z
	double tmp = LERP(tmp1234, tmp5678, fractz);
	return tmp;
  #endif
}

vec3 GridData::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0] - theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1] - theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2] - theCellSize*0.5), mMax[2]);
   return out;
}

GridDataX::GridDataX() : GridData()
{
}

GridDataX::~GridDataX()
{
}

void GridDataX::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*(theDim[0]+1);
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*theDim[2];
   mData.resize((theDim[0]+1)*theDim[1]*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataX::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (i < 0 || i > theDim[0]) return dflt;

   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*(theDim[0]+1);
   int stack = j*(theDim[0]+1)*theDim[2];
   return mData[stack + row + col];
}

const double GridDataX::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (i < 0 || i > theDim[0]) return dflt;

   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*(theDim[0]+1);
   int stack = j*(theDim[0]+1)*theDim[2];
   return mData[stack + row + col];
}

vec3 GridDataX::worldToSelf(const vec3& pt) const
{   
   vec3 out;
   out[0] = min(max(0.0, pt[0]), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataY::GridDataY() : GridData()
{
}

GridDataY::~GridDataY()
{
}

void GridDataY::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*(theDim[1]+1);
   mMax[2] = theCellSize*theDim[2];
   mData.resize(theDim[0]*(theDim[1]+1)*theDim[2], false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataY::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (j < 0 || j > theDim[1]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return mData[stack + row + col];
}

const double GridDataY::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (j < 0 || j > theDim[1]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (k < 0) k = 0;
   if (k > theDim[2]-1) k = theDim[2]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*theDim[2];
   return mData[stack + row + col];
}

vec3 GridDataY::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]), mMax[1]);
   out[2] = min(max(0.0, pt[2]-theCellSize*0.5), mMax[2]);
   return out;
}

GridDataZ::GridDataZ() : GridData()
{
}

GridDataZ::~GridDataZ()
{
}

void GridDataZ::initialize(double dfltValue)
{
   GridData::initialize(dfltValue);
   mMax[0] = theCellSize*theDim[0];
   mMax[1] = theCellSize*theDim[1];
   mMax[2] = theCellSize*(theDim[2]+1);
   mData.resize(theDim[0]*theDim[1]*(theDim[2]+1), false);
   std::fill(mData.begin(), mData.end(), mDfltValue);
}

double& GridDataZ::operator()(int i, int j, int k)
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (k < 0 || k > theDim[2]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*(theDim[2]+1);

   return mData[stack + row + col];
}

const double GridDataZ::operator()(int i, int j, int k) const
{
   static double dflt = 0;
   dflt = mDfltValue;  // Protect against setting the default value

   if (k < 0 || k > theDim[2]) return dflt;

   if (i < 0) i = 0;
   if (i > theDim[0]-1) i = theDim[0]-1;
   if (j < 0) j = 0;
   if (j > theDim[1]-1) j = theDim[1]-1;

   int col = i;
   int row = k*theDim[0];
   int stack = j*theDim[0]*(theDim[2]+1);

   return mData[stack + row + col];
}

vec3 GridDataZ::worldToSelf(const vec3& pt) const
{
   vec3 out;
   out[0] = min(max(0.0, pt[0]-theCellSize*0.5), mMax[0]);
   out[1] = min(max(0.0, pt[1]-theCellSize*0.5), mMax[1]);
   out[2] = min(max(0.0, pt[2]), mMax[2]);
   return out;
}
