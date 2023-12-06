
//  Copyright (c) 2023, Andre Caceres Carrilho
//
//  Mathlib.h common linear algebra algorithms
//  this code needs at least C++17

#ifndef _Mathlib_h 
#define _Mathlib_h 1

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

namespace Math { 

template<typename Real>
inline Real degToRad(const Real& x)
{ return x * Real(3.1415926535897932384626) / Real(180); }

template<typename Real>
inline Real radToDeg(const Real& x)
{ return x * Real(180) / Real(3.1415926535897932384626); }

template<typename Real>
constexpr Real eps()  // Machine epsilon
{ return std::numeric_limits<Real>::epsilon() * Real(100); }

// According to IEEE-754: eps = 2^(digits-1)
//  machine epsilon for float  2^-23 ~ 1.19 e-07  
//  machine epsilon for double 2^-52 ~ 2.22 e-16   
//
// Please note that MSVC treats long double == double (64 bit), hence 
// the same epsilon is hard-coded here. GCC uses a 80 bit long double.
template<> constexpr float       eps<float>()       { return 1e-5f; }
template<> constexpr double      eps<double>()      { return 1e-11; }
template<> constexpr long double eps<long double>() { return 1e-11; }

#include "DataLayout.h"
#include "Matrix.h"
#include "Operations.h"
#include "Cholesky.h"
#include "QR.h"
#include "LU.h"
#include "SVD.h"
#include "Eigen.h"
#include "EigenSymmetric.h"
#include "EigenRegular.h"
#include "Functions.h"

typedef MatrixGeneric<DenseRowMajor<double>> Matrix;

typedef MatrixGeneric<DenseColMajorFixed<double,2,2>> Matrix2d;
typedef MatrixGeneric<DenseColMajorFixed<double,2,1>> Vector2d;
typedef MatrixGeneric<DenseColMajorFixed<float,2,2>> Matrix2f;
typedef MatrixGeneric<DenseColMajorFixed<float,2,1>> Vector2f;

typedef MatrixGeneric<DenseColMajorFixed<double,3,3>> Matrix3d;
typedef MatrixGeneric<DenseColMajorFixed<double,3,1>> Vector3d;
typedef MatrixGeneric<DenseColMajorFixed<float,3,3>> Matrix3f;
typedef MatrixGeneric<DenseColMajorFixed<float,3,1>> Vector3f;

typedef MatrixGeneric<DenseColMajorFixed<double,4,4>> Matrix4d;
typedef MatrixGeneric<DenseColMajorFixed<double,4,1>> Vector4d;
typedef MatrixGeneric<DenseColMajorFixed<float,4,4>> Matrix4f;
typedef MatrixGeneric<DenseColMajorFixed<float,4,1>> Vector4f;

} // Math

#endif 
