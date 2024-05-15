
/// \author Andre Caceres Carrilho
/// \namespace math  
/// \file math3d.cc

#include <cmath>
#include <algorithm>
#include "math3d.h"
namespace math{

/////////////////////////////////////////////////////////////////
// float 
/////////////////////////////////////////////////////////////////
Vector4f::Vector4f(float X, float Y, float Z, float W)
#ifdef _SIMD_x86
{ r = _mm_setr_ps(X, Y, Z, W); }
#else 
    : x(X)
    , y(Y)
    , z(Z)
    , w(W)
{}
#endif 

Vector4f::Vector4f(const Vector4f &a)
#ifdef _SIMD_x86
{ r = a.r; }
#else 
    : x(a.x)
    , y(a.y)
    , z(a.z)
    , w(a.w)
{}
#endif 

Vector4f& Vector4f::operator= (const Vector4f& a) {
#ifdef _SIMD_x86
    r = a.r;
#else 
    x = a.x;  
    y = a.y;  
    z = a.z;  
    w = a.w;
#endif
    return *this;
}

Vector4f& Vector4f::operator= (float s) {
#ifdef _SIMD_x86
    r = _mm_set1_ps(s);
#else 
    x = s; 
    y = s; 
    z = s; 
    w = s;
#endif 
    return *this;
}

Vector4f& Vector4f::operator+= (const Vector4f& a) {
#ifdef _SIMD_x86
    r = _mm_add_ps(r, a.r);
#else
    x += a.x;  
    y += a.y;  
    z += a.z;  
    w += a.w;
#endif
    return *this;
}

Vector4f& Vector4f::operator-= (const Vector4f& a) {
#ifdef _SIMD_x86
    r = _mm_sub_ps(r, a.r);
#else
    x -= a.x;  
    y -= a.y;  
    z -= a.z;  
    w -= a.w;
#endif
    return *this;
}

/// \brief Element-wise product 
/// \param a Vector with elements to be multiplied
/// \return this
Vector4f& Vector4f::operator*= (const Vector4f& a) {
#ifdef _SIMD_x86
    r = _mm_mul_ps(r, a.r);
#else
    x *= a.x;  
    y *= a.y;  
    z *= a.z;  
    w *= a.w;
#endif
    return *this;
}

Vector4f& Vector4f::operator*= (float s) {
#ifdef _SIMD_x86
    r = _mm_mul_ps(r, _mm_set1_ps(s));
#else 
    x *= s; 
    y *= s; 
    z *= s; 
    w *= s;
#endif
    return *this;
}

float Vector4f::norm() const {
    return std::sqrt(x*x + y*y + z*z);
}

void Vector4f::normalize(){
    (*this) *= 1.0f / norm();
}

void Vector4f::homogenize() {
#ifdef _SIMD_x86
    r = _mm_div_ps(r, _mm_set1_ps(w));
#else
    x /= w;  
    y /= w;  
    z /= w;  
    w /= w;
#endif 
}

Vector4f& Vector4f::max(const Vector4f &a) {
#ifdef _SIMD_x86
    r = _mm_max_ps(r, a.r);
#else 
    x = std::max(x, a.x);
    y = std::max(y, a.y);
    z = std::max(z, a.z);
    w = std::max(w, a.w);
#endif 
    return *this;
}

Vector4f& Vector4f::min(const Vector4f &a) {
#ifdef _SIMD_x86
    r = _mm_min_ps(r, a.r);
#else 
    x = std::min(x, a.x);
    y = std::min(y, a.y);
    z = std::min(z, a.z);
    w = std::min(w, a.w);
#endif 
    return *this;
}

bool Vector4f::operator > (const Vector4f &a) const {
#ifdef _SIMD_x86
    return _mm_movemask_ps(_mm_cmp_ps(r, a.r, _CMP_GT_OQ)) == 0xf;
#else 
    return x > a.x 
        && y > a.y 
        && z > a.z 
        && w > a.w;
#endif 
}

bool Vector4f::operator < (const Vector4f &a) const {
#ifdef _SIMD_x86
    return _mm_movemask_ps(_mm_cmp_ps(r, a.r, _CMP_LT_OQ)) == 0xf;
#else 
    return x < a.x 
        && y < a.y 
        && z < a.z 
        && w < a.w;
#endif 
}

bool Vector4f::operator >= (const Vector4f &a) const {
#ifdef _SIMD_x86
    return _mm_movemask_ps(_mm_cmp_ps(r, a.r, _CMP_GE_OQ)) == 0xf;
#else 
    return x >= a.x 
        && y >= a.y 
        && z >= a.z 
        && w >= a.w;
#endif 
}

bool Vector4f::operator <= (const Vector4f &a) const {
#ifdef _SIMD_x86
    return _mm_movemask_ps(_mm_cmp_ps(r, a.r, _CMP_LE_OQ)) == 0xf;
#else 
    return x <= a.x 
        && y <= a.y 
        && z <= a.z 
        && w <= a.w;
#endif 
}


/////////////////////////////////////////////////////////////////
// double 
/////////////////////////////////////////////////////////////////
Vector4d::Vector4d(double X, double Y, double Z, double W)
#ifdef _SIMD_x86
{ r = _mm256_setr_pd(X, Y, Z, W); }
#else 
    : x(X)
    , y(Y)
    , z(Z)
    , w(W)
{}
#endif 

Vector4d::Vector4d(const Vector4d &a)
#ifdef _SIMD_x86
{ r = a.r; }
#else 
    : x(a.x)
    , y(a.y)
    , z(a.z)
    , w(a.w)
{}
#endif 

Vector4d& Vector4d::operator= (const Vector4d& a) {
#ifdef _SIMD_x86
    r = a.r;
#else 
    x = a.x;  
    y = a.y;  
    z = a.z;  
    w = a.w;
#endif
    return *this;
}

Vector4d& Vector4d::operator= (double s) {
#ifdef _SIMD_x86
    r = _mm256_set1_pd(s);
#else 
    x = s; 
    y = s; 
    z = s; 
    w = s;
#endif 
    return *this;
}

Vector4d& Vector4d::operator+= (const Vector4d& a) {
#ifdef _SIMD_x86
    r = _mm256_add_pd(r, a.r);
#else
    x += a.x;  
    y += a.y;  
    z += a.z;  
    w += a.w;
#endif
    return *this;
}

Vector4d& Vector4d::operator-= (const Vector4d& a) {
#ifdef _SIMD_x86
    r = _mm256_sub_pd(r, a.r);
#else
    x -= a.x;  
    y -= a.y;  
    z -= a.z;  
    w -= a.w;
#endif
    return *this;
}

/// \brief Element-wise product 
/// \param a Vector with elements to be multiplied
/// \return this
Vector4d& Vector4d::operator*= (const Vector4d& a) {
#ifdef _SIMD_x86
    r = _mm256_mul_pd(r, a.r);
#else
    x *= a.x;  
    y *= a.y;  
    z *= a.z;  
    w *= a.w;
#endif
    return *this;
}

Vector4d& Vector4d::operator*= (double s) {
#ifdef _SIMD_x86
    r = _mm256_mul_pd(r, _mm256_set1_pd(s));
#else 
    x *= s; 
    y *= s; 
    z *= s; 
    w *= s;
#endif
    return *this;
}

/// \brief Divide all elements by w: \n
/// \f$ x/=w, y/=w, z/=w, w/=w \f$
/// \attention Doesn't check if w is not null
void Vector4d::homogenize() {
#ifdef _SIMD_x86
    r = _mm256_div_pd(r, _mm256_set1_pd(w));
#else
    x /= w;  
    y /= w;  
    z /= w;  
    w /= w;
#endif 
}

Vector4d& Vector4d::max(const Vector4d &a) {
#ifdef _SIMD_x86
    r = _mm256_max_pd(r, a.r);
#else 
    x = std::max(x, a.x);
    y = std::max(y, a.y);
    z = std::max(z, a.z);
    w = std::max(w, a.w);
#endif 
    return *this;
}

Vector4d& Vector4d::min(const Vector4d &a) {
#ifdef _SIMD_x86
    r = _mm256_min_pd(r, a.r);
#else 
    x = std::min(x, a.x);
    y = std::min(y, a.y);
    z = std::min(z, a.z);
    w = std::min(w, a.w);
#endif 
    return *this;
}

bool Vector4d::operator > (const Vector4d &a) const {
#ifdef _SIMD_x86
    return _mm256_movemask_pd(_mm256_cmp_pd(r, a.r, _CMP_GT_OQ)) == 0xf;
#else 
    return x > a.x 
        && y > a.y 
        && z > a.z 
        && w > a.w;
#endif 
}

bool Vector4d::operator < (const Vector4d &a) const {
#ifdef _SIMD_x86
    return _mm256_movemask_pd(_mm256_cmp_pd(r, a.r, _CMP_LT_OQ)) == 0xf;
#else 
    return x < a.x 
        && y < a.y 
        && z < a.z 
        && w < a.w;
#endif 
}

bool Vector4d::operator >= (const Vector4d &a) const {
#ifdef _SIMD_x86
    return _mm256_movemask_pd(_mm256_cmp_pd(r, a.r, _CMP_GE_OQ)) == 0xf;
#else 
    return x >= a.x 
        && y >= a.y 
        && z >= a.z 
        && w >= a.w;
#endif 
}

bool Vector4d::operator <= (const Vector4d &a) const {
#ifdef _SIMD_x86
    return _mm256_movemask_pd(_mm256_cmp_pd(r, a.r, _CMP_LE_OQ)) == 0xf;
#else 
    return x <= a.x 
        && y <= a.y 
        && z <= a.z 
        && w <= a.w;
#endif 
}

/// \brief Computes the 3D norm of the Vector
/// \return \f$ sqrt(x*x + y*y + z*z) \f$
/// \attention only (X,Y,Z) are used in the formula!
double Vector4d::norm() const {
    return std::sqrt(x*x + y*y + z*z);
}

/// \brief Sets 'this' Vector length to one, i.e. \f$ sqrt(x*x + y*y + z*z) = 1 \f$
/// \note W is garbage after
/// \attention Doesn't check if the Vector norm is not null 
void Vector4d::normalize(){
    (*this) *= 1.0 / norm();
}

/// \brief Squared Euclidian distance between \p a and \p b 
/// \param a  3D point
/// \param b  3D point
/// \return \f$ dx^2 + dy^2 + dz^2 \f$
/// \attention only (X,Y,Z) are used in the formula!
double dist2(const Vector4d &a, const Vector4d &b){
    Vector4d v(a);
    v -= b;
    return dot(v, v);
}
float dist2(const Vector4f &a, const Vector4f &b){
    Vector4f v(a);
    v -= b;
    return dot(v, v);
}

/// \brief Euclidian distance between \p a and \p b 
/// \param a  3D point
/// \param b  3D point
/// \return \f$ sqrt( dx^2 + dy^2 + dz^2 ) \f$
/// \attention only (X,Y,Z) are used in the formula!
double dist(const Vector4d &a, const Vector4d &b){
    return std::sqrt(dist(a,b));
}
float dist(const Vector4f &a, const Vector4f &b){
    return std::sqrt(dist(a,b));
}

/// \brief Dot Product between Vector \p a and \p b
/// \param a 3D Vector
/// \param b 3D Vector
/// \return ax*bx + ay*by + az*bz
/// \attention only (X,Y,Z) are used in the formula!
double dot(const Vector4d &a, const Vector4d &b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
float dot(const Vector4f &a, const Vector4f &b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/// \brief Cross Product between Vector \p a and \p b
/// \note \f$ c = a * b \f$
/// \attention only (X,Y,Z) are used in the formula!
void cross(Vector4d &c, const Vector4d &a, const Vector4d &b){
#ifdef _SIMD_x86
    c.r = _mm256_permute4x64_pd(
        _mm256_sub_pd(
			_mm256_mul_pd(a.r, _mm256_permute4x64_pd(b.r, _MM_SHUFFLE(3, 0, 2, 1))),
			_mm256_mul_pd(b.r, _mm256_permute4x64_pd(a.r, _MM_SHUFFLE(3, 0, 2, 1)))
        ), _MM_SHUFFLE(3, 0, 2, 1)
    );
#else 
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
#endif 
}
void cross(Vector4f &c, const Vector4f &a, const Vector4f &b){
#ifdef _SIMD_x86    
    c.r = _mm_permute_ps(
        _mm_sub_ps(
			_mm_mul_ps(a.r, _mm_permute_ps(b.r, _MM_SHUFFLE(3, 0, 2, 1))),
			_mm_mul_ps(b.r, _mm_permute_ps(a.r, _MM_SHUFFLE(3, 0, 2, 1)))
        ), _MM_SHUFFLE(3, 0, 2, 1)
    );
#else 
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
#endif 
}

void corr(Matrix4f &C, const Vector4f &a, const Vector4f &b){
#ifdef _SIMD_x86
    C.r0 = _mm_mul_ps(_mm_set1_ps(a.x), b.r);
    C.r1 = _mm_mul_ps(_mm_set1_ps(a.y), b.r);
    C.r2 = _mm_mul_ps(_mm_set1_ps(a.z), b.r);
    C.r3 = _mm_mul_ps(_mm_set1_ps(a.w), b.r);
#else 
    C.m[ 0]=a.x*b.x; 
    C.m[ 1]=a.x*b.y;
    C.m[ 2]=a.x*b.z;
    C.m[ 3]=a.x*b.w;

    C.m[ 4]=a.y*b.x; 
    C.m[ 5]=a.y*b.y;
    C.m[ 6]=a.y*b.z;
    C.m[ 7]=a.y*b.w;
    
    C.m[ 8]=a.z*b.x; 
    C.m[ 9]=a.z*b.y;
    C.m[10]=a.z*b.z;
    C.m[11]=a.z*b.w;

    C.m[12]=a.w*b.x; 
    C.m[13]=a.w*b.y;
    C.m[14]=a.w*b.z;
    C.m[15]=a.w*b.w;
#endif 
}

void corr(Matrix4f &C, const Vector4f &a){
#ifdef _SIMD_x86
    C.r0 = _mm_mul_ps(_mm_set1_ps(a.x), a.r);
    C.r1 = _mm_mul_ps(_mm_set1_ps(a.y), a.r);
    C.r2 = _mm_mul_ps(_mm_set1_ps(a.z), a.r);
    C.r3 = _mm_mul_ps(_mm_set1_ps(a.w), a.r);
#else 
#if 1
    C.m[ 0]=a.x*a.x; 
    C.m[ 1]=a.x*a.y;
    C.m[ 2]=a.x*a.z;
    C.m[ 3]=a.x*a.w;
  
    C.m[ 4]=C.m[ 1]; //  C.m[ 4]=a.y*a.x; 
    C.m[ 5]=a.y*a.y;
    C.m[ 6]=a.y*a.z;
    C.m[ 7]=a.y*a.w;
    
    C.m[ 8]=C.m[ 2]; //  C.m[ 8]=a.z*a.x; 
    C.m[ 9]=C.m[ 6]; //  C.m[ 9]=a.z*a.y;
    C.m[10]=a.z*a.z;
    C.m[11]=a.z*a.w;

    C.m[12]=C.m[ 3]; //  C.m[12]=a.w*a.x; 
    C.m[13]=C.m[ 7]; //  C.m[13]=a.w*a.y;
    C.m[14]=C.m[11]; //  C.m[14]=a.w*a.z;
    C.m[15]=a.w*a.w;
#else 
    C.m[ 0]=a.x*a.x; 
    C.m[ 1]=a.x*a.y;
    C.m[ 2]=a.x*a.z;
    C.m[ 3]=a.x*a.w;
  
    C.m[ 4]=a.y*a.x; 
    C.m[ 5]=a.y*a.y;
    C.m[ 6]=a.y*a.z;
    C.m[ 7]=a.y*a.w;
    
    C.m[ 8]=a.z*a.x; 
    C.m[ 9]=a.z*a.y;
    C.m[10]=a.z*a.z;
    C.m[11]=a.z*a.w;

    C.m[12]=a.w*a.x;  
    C.m[13]=a.w*a.y; 
    C.m[14]=a.w*a.z; 
    C.m[15]=a.w*a.w;
#endif 
#endif 
}

void corr(Matrix4d &C, const Vector4d &a){
#ifdef _SIMD_x86
    C.r0 = _mm256_mul_pd(_mm256_set1_pd(a.x), a.r);
    C.r1 = _mm256_mul_pd(_mm256_set1_pd(a.y), a.r);
    C.r2 = _mm256_mul_pd(_mm256_set1_pd(a.z), a.r);
    C.r3 = _mm256_mul_pd(_mm256_set1_pd(a.w), a.r);
#else 
#if 1
    C.m[ 0]=a.x*a.x; 
    C.m[ 1]=a.x*a.y;
    C.m[ 2]=a.x*a.z;
    C.m[ 3]=a.x*a.w;
  
    C.m[ 4]=C.m[ 1]; //  C.m[ 4]=a.y*a.x; 
    C.m[ 5]=a.y*a.y;
    C.m[ 6]=a.y*a.z;
    C.m[ 7]=a.y*a.w;
    
    C.m[ 8]=C.m[ 2]; //  C.m[ 8]=a.z*a.x; 
    C.m[ 9]=C.m[ 6]; //  C.m[ 9]=a.z*a.y;
    C.m[10]=a.z*a.z;
    C.m[11]=a.z*a.w;

    C.m[12]=C.m[ 3]; //  C.m[12]=a.w*a.x; 
    C.m[13]=C.m[ 7]; //  C.m[13]=a.w*a.y;
    C.m[14]=C.m[11]; //  C.m[14]=a.w*a.z;
    C.m[15]=a.w*a.w;
#else 
    C.m[ 0]=a.x*a.x; 
    C.m[ 1]=a.x*a.y;
    C.m[ 2]=a.x*a.z;
    C.m[ 3]=a.x*a.w;
  
    C.m[ 4]=a.y*a.x; 
    C.m[ 5]=a.y*a.y;
    C.m[ 6]=a.y*a.z;
    C.m[ 7]=a.y*a.w;
    
    C.m[ 8]=a.z*a.x; 
    C.m[ 9]=a.z*a.y;
    C.m[10]=a.z*a.z;
    C.m[11]=a.z*a.w;

    C.m[12]=a.w*a.x;  
    C.m[13]=a.w*a.y; 
    C.m[14]=a.w*a.z; 
    C.m[15]=a.w*a.w;
#endif 
#endif 
}

void corr(Matrix4d &C, const Vector4d &a, const Vector4d &b){
#ifdef _SIMD_x86
    C.r0 = _mm256_mul_pd(_mm256_set1_pd(a.x), b.r);
    C.r1 = _mm256_mul_pd(_mm256_set1_pd(a.y), b.r);
    C.r2 = _mm256_mul_pd(_mm256_set1_pd(a.z), b.r);
    C.r3 = _mm256_mul_pd(_mm256_set1_pd(a.w), b.r);
#else 
    C.m[ 0]=a.x*b.x; 
    C.m[ 1]=a.x*b.y;
    C.m[ 2]=a.x*b.z;
    C.m[ 3]=a.x*b.w;

    C.m[ 4]=a.y*b.x; 
    C.m[ 5]=a.y*b.y;
    C.m[ 6]=a.y*b.z;
    C.m[ 7]=a.y*b.w;
    
    C.m[ 8]=a.z*b.x; 
    C.m[ 9]=a.z*b.y;
    C.m[10]=a.z*b.z;
    C.m[11]=a.z*b.w;

    C.m[12]=a.w*b.x; 
    C.m[13]=a.w*b.y;
    C.m[14]=a.w*b.z;
    C.m[15]=a.w*b.w;
#endif 
}

/////////////////////////////////////////////////////////////////
// float 
/////////////////////////////////////////////////////////////////

///       |  0   1   2   3 |
///  A =  |  4   5   6   7 |
///       |  8   9  10  11 |
///       | 12  13  14  15 |


Matrix4f::Matrix4f(){
#ifdef _SIMD_x86
    r0 = _mm_setzero_ps();
    r1 = _mm_setzero_ps();
    r2 = _mm_setzero_ps();
    r3 = _mm_setzero_ps();
#else 
    m[ 0]=0.0f; m[ 1]=0.0f; m[ 2]=0.0f; m[ 3]=0.0f;
    m[ 4]=0.0f; m[ 5]=0.0f; m[ 6]=0.0f; m[ 7]=0.0f;
    m[ 8]=0.0f; m[ 9]=0.0f; m[10]=0.0f; m[11]=0.0f;
    m[12]=0.0f; m[13]=0.0f; m[14]=0.0f; m[15]=0.0f;
#endif 
}

Matrix4f::Matrix4f(const Matrix4f &A){
#ifdef _SIMD_x86
    r0 = A.r0;
    r1 = A.r1;
    r2 = A.r2;
    r3 = A.r3;
#else 
    m[ 0]=A.m[ 0]; m[ 1]=A.m[ 1]; m[ 2]=A.m[ 2]; m[ 3]=A.m[ 3];
    m[ 4]=A.m[ 4]; m[ 5]=A.m[ 5]; m[ 6]=A.m[ 6]; m[ 7]=A.m[ 7];
    m[ 8]=A.m[ 8]; m[ 9]=A.m[ 9]; m[10]=A.m[10]; m[11]=A.m[11];
    m[12]=A.m[12]; m[13]=A.m[13]; m[14]=A.m[14]; m[15]=A.m[15];
#endif 
}

Matrix4f &Matrix4f::operator = (const Matrix4f &A){
#ifdef _SIMD_x86
    r0 = A.r0;
    r1 = A.r1;
    r2 = A.r2;
    r3 = A.r3;
#else 
    m[ 0]=A.m[ 0]; m[ 1]=A.m[ 1]; m[ 2]=A.m[ 2]; m[ 3]=A.m[ 3];
    m[ 4]=A.m[ 4]; m[ 5]=A.m[ 5]; m[ 6]=A.m[ 6]; m[ 7]=A.m[ 7];
    m[ 8]=A.m[ 8]; m[ 9]=A.m[ 9]; m[10]=A.m[10]; m[11]=A.m[11];
    m[12]=A.m[12]; m[13]=A.m[13]; m[14]=A.m[14]; m[15]=A.m[15];
#endif 
    return *this;
}

Matrix4f& Matrix4f::operator = (float s){
#ifdef _SIMD_x86
    r0 = _mm_set1_ps(s);
    r1 = _mm_set1_ps(s);
    r2 = _mm_set1_ps(s);
    r3 = _mm_set1_ps(s);
#else 
    m[ 0]=s; m[ 1]=s; m[ 2]=s; m[ 3]=s;
    m[ 4]=s; m[ 5]=s; m[ 6]=s; m[ 7]=s;
    m[ 8]=s; m[ 9]=s; m[10]=s; m[11]=s;
    m[12]=s; m[13]=s; m[14]=s; m[15]=s;
#endif 
    return *this;
}

Matrix4f& Matrix4f::operator += (const Matrix4f &A){
#ifdef _SIMD_x86
    r0 = _mm_add_ps(r0, A.r0);
    r1 = _mm_add_ps(r1, A.r1);
    r2 = _mm_add_ps(r2, A.r2);
    r3 = _mm_add_ps(r3, A.r3);
#else 
    m[ 0]+=A.m[ 0]; m[ 1]+=A.m[ 1]; m[ 2]+=A.m[ 2]; m[ 3]+=A.m[ 3];
    m[ 4]+=A.m[ 4]; m[ 5]+=A.m[ 5]; m[ 6]+=A.m[ 6]; m[ 7]+=A.m[ 7];
    m[ 8]+=A.m[ 8]; m[ 9]+=A.m[ 9]; m[10]+=A.m[10]; m[11]+=A.m[11];
    m[12]+=A.m[12]; m[13]+=A.m[13]; m[14]+=A.m[14]; m[15]+=A.m[15];
#endif 
    return *this;
}

Matrix4f& Matrix4f::operator -= (const Matrix4f &A){
#ifdef _SIMD_x86
    r0 = _mm_sub_ps(r0, A.r0);
    r1 = _mm_sub_ps(r1, A.r1);
    r2 = _mm_sub_ps(r2, A.r2);
    r3 = _mm_sub_ps(r3, A.r3);
#else 
    m[ 0]-=A.m[ 0]; m[ 1]-=A.m[ 1]; m[ 2]-=A.m[ 2]; m[ 3]-=A.m[ 3];
    m[ 4]-=A.m[ 4]; m[ 5]-=A.m[ 5]; m[ 6]-=A.m[ 6]; m[ 7]-=A.m[ 7];
    m[ 8]-=A.m[ 8]; m[ 9]-=A.m[ 9]; m[10]-=A.m[10]; m[11]-=A.m[11];
    m[12]-=A.m[12]; m[13]-=A.m[13]; m[14]-=A.m[14]; m[15]-=A.m[15];
#endif 
    return *this;
}

Matrix4f& Matrix4f::operator *= (float s){
#ifdef _SIMD_x86
    __m128 vs = _mm_set1_ps(s);
    r0 = _mm_mul_ps(r0, vs);
    r1 = _mm_mul_ps(r1, vs);
    r2 = _mm_mul_ps(r2, vs);
    r3 = _mm_mul_ps(r3, vs);
#else 
    m[ 0]*=s; m[ 1]*=s; m[ 2]*=s; m[ 3]*=s;
    m[ 4]*=s; m[ 5]*=s; m[ 6]*=s; m[ 7]*=s;
    m[ 8]*=s; m[ 9]*=s; m[10]*=s; m[11]*=s;
    m[12]*=s; m[13]*=s; m[14]*=s; m[15]*=s;
#endif 
    return *this;
}

Matrix4f& Matrix4f::operator *= (const Matrix4f &A){
#ifdef _SIMD_x86
    __m128 row;
    Vector4f t;

    t.r = r0;
    row = _mm_mul_ps  (A.r0, _mm_set1_ps(t.x));
    row = _mm_fmadd_ps(A.r1, _mm_set1_ps(t.y), row);
    row = _mm_fmadd_ps(A.r2, _mm_set1_ps(t.z), row);
    row = _mm_fmadd_ps(A.r3, _mm_set1_ps(t.w), row);
    r0 = row;

    t.r = r1;
    row = _mm_mul_ps  (A.r0, _mm_set1_ps(t.x));
    row = _mm_fmadd_ps(A.r1, _mm_set1_ps(t.y), row);
    row = _mm_fmadd_ps(A.r2, _mm_set1_ps(t.z), row);
    row = _mm_fmadd_ps(A.r3, _mm_set1_ps(t.w), row);
    r1 = row;

    t.r = r2;
    row = _mm_mul_ps  (A.r0, _mm_set1_ps(t.x));
    row = _mm_fmadd_ps(A.r1, _mm_set1_ps(t.y), row);
    row = _mm_fmadd_ps(A.r2, _mm_set1_ps(t.z), row);
    row = _mm_fmadd_ps(A.r3, _mm_set1_ps(t.w), row);
    r2 = row;
    
    t.r = r3;
    row = _mm_mul_ps  (A.r0, _mm_set1_ps(t.x));
    row = _mm_fmadd_ps(A.r1, _mm_set1_ps(t.y), row);
    row = _mm_fmadd_ps(A.r2, _mm_set1_ps(t.z), row);
    row = _mm_fmadd_ps(A.r3, _mm_set1_ps(t.w), row);
    r3 = row;
#else 
    float t[4];

    t[0]=m[ 0]; t[1]=m[ 1]; t[2]=m[ 2]; t[3]=m[ 3];
    m[ 0]=t[0]*A.m[ 0] + t[1]*A.m[ 4] + t[2]*A.m[ 8] + t[3]*A.m[12];
    m[ 1]=t[0]*A.m[ 1] + t[1]*A.m[ 5] + t[2]*A.m[ 9] + t[3]*A.m[13];
    m[ 2]=t[0]*A.m[ 2] + t[1]*A.m[ 6] + t[2]*A.m[10] + t[3]*A.m[14];
    m[ 3]=t[0]*A.m[ 3] + t[1]*A.m[ 7] + t[2]*A.m[11] + t[3]*A.m[15];

    t[0]=m[ 4]; t[1]=m[ 5]; t[2]=m[ 6]; t[3]=m[ 7];
    m[ 4]=t[0]*A.m[ 0] + t[1]*A.m[ 4] + t[2]*A.m[ 8] + t[3]*A.m[12];
    m[ 5]=t[0]*A.m[ 1] + t[1]*A.m[ 5] + t[2]*A.m[ 9] + t[3]*A.m[13];
    m[ 6]=t[0]*A.m[ 2] + t[1]*A.m[ 6] + t[2]*A.m[10] + t[3]*A.m[14];
    m[ 7]=t[0]*A.m[ 3] + t[1]*A.m[ 7] + t[2]*A.m[11] + t[3]*A.m[15];

    t[0]=m[ 8]; t[1]=m[ 9]; t[2]=m[10]; t[3]=m[11];
    m[ 8]=t[0]*A.m[ 0] + t[1]*A.m[ 4] + t[2]*A.m[ 8] + t[3]*A.m[12];
    m[ 9]=t[0]*A.m[ 1] + t[1]*A.m[ 5] + t[2]*A.m[ 9] + t[3]*A.m[13];
    m[10]=t[0]*A.m[ 2] + t[1]*A.m[ 6] + t[2]*A.m[10] + t[3]*A.m[14];
    m[11]=t[0]*A.m[ 3] + t[1]*A.m[ 7] + t[2]*A.m[11] + t[3]*A.m[15];

    t[0]=m[12]; t[1]=m[13]; t[2]=m[14]; t[3]=m[15];
    m[12]=t[0]*A.m[ 0] + t[1]*A.m[ 4] + t[2]*A.m[ 8] + t[3]*A.m[12];
    m[13]=t[0]*A.m[ 1] + t[1]*A.m[ 5] + t[2]*A.m[ 9] + t[3]*A.m[13];
    m[14]=t[0]*A.m[ 2] + t[1]*A.m[ 6] + t[2]*A.m[10] + t[3]*A.m[14];
    m[15]=t[0]*A.m[ 3] + t[1]*A.m[ 7] + t[2]*A.m[11] + t[3]*A.m[15];
#endif
    return *this;
}

void Matrix4f::transform(Vector4f &a, const Vector4f &b){
#ifdef _SIMD_x86
    __m128 k;
    __m128 col0 = _mm_setr_ps(m[0], m[4], m[ 8], m[12]);
    __m128 col1 = _mm_setr_ps(m[1], m[5], m[ 9], m[13]);
    __m128 col2 = _mm_setr_ps(m[2], m[6], m[10], m[14]);
    __m128 col3 = _mm_setr_ps(m[3], m[7], m[11], m[15]);

    k = _mm_mul_ps  (col0, _mm_set1_ps(b.x));
    k = _mm_fmadd_ps(col1, _mm_set1_ps(b.y), k);
    k = _mm_fmadd_ps(col2, _mm_set1_ps(b.z), k);
    k = _mm_fmadd_ps(col3, _mm_set1_ps(b.w), k);

    a.r = k;
#else
    a.v[0] = m[ 0]*b.v[0] + m[ 1]*b.v[1] + m[ 2]*b.v[2] + m[ 3]*b.v[3];
    a.v[1] = m[ 4]*b.v[0] + m[ 5]*b.v[1] + m[ 6]*b.v[2] + m[ 7]*b.v[3];
    a.v[2] = m[ 8]*b.v[0] + m[ 9]*b.v[1] + m[10]*b.v[2] + m[11]*b.v[3];
    a.v[3] = m[12]*b.v[0] + m[13]*b.v[1] + m[14]*b.v[2] + m[15]*b.v[3];
#endif
}

void Matrix4f::invert()
{
    Vector4f tu;
    Vector4f tv;

    /// 2x2 sub-derminants
    Vector4f d21;
    tu = Vector4f(m[ 0], m[ 0], m[ 0], m[ 1]);
    tv = Vector4f(m[ 5], m[ 6], m[ 7], m[ 6]);
    tu *= tv;
    d21 = tu;
    tu = Vector4f(m[ 1], m[ 2], m[ 3], m[ 2]);
    tv = Vector4f(m[ 4], m[ 4], m[ 4], m[ 5]);
    tu *= tv;
    d21 -= tu;
    
    Vector4f d22;
    tu = Vector4f(m[ 1], m[ 2], m[ 0], m[ 0]);
    tv = Vector4f(m[ 7], m[ 7], m[13], m[14]);
    tu *= tv;
    d22 = tu;
    tu = Vector4f(m[ 3], m[ 3], m[ 1], m[ 2]);
    tv = Vector4f(m[ 5], m[ 6], m[12], m[12]);
    tu *= tv;
    d22 -= tu;
    
    Vector4f d23;
    tu = Vector4f(m[ 0], m[ 1], m[ 1], m[ 2]);
    tv = Vector4f(m[15], m[14], m[15], m[15]);
    tu *= tv;
    d23 = tu;
    tu = Vector4f(m[ 3], m[ 2], m[ 3], m[ 3]);
    tv = Vector4f(m[12], m[13], m[13], m[14]);
    tu *= tv;
    d23 -= tu;
    
    Vector4f d24;
    tu = Vector4f(m[ 4], m[ 4], m[ 4], m[ 5]);
    tv = Vector4f(m[13], m[14], m[15], m[14]);
    tu *= tv;
    d24 = tu;
    tu = Vector4f(m[ 5], m[ 6], m[ 7], m[ 6]);
    tv = Vector4f(m[12], m[12], m[12], m[13]);
    tu *= tv;
    d24 -= tu;
    
    Vector4f d25;
    tu = Vector4f(m[ 5], m[ 6], 0, 0);
    tv = Vector4f(m[15], m[15], 0, 0);
    tu *= tv;
    d25 = tu;
    tu = Vector4f(m[ 7], m[ 7], 0, 0);
    tv = Vector4f(m[13], m[14], 0, 0);
    tu *= tv;
    d25 -= tu;

    /// temp 
    Vector4f r2 = v[2];
    Vector4f r3 = v[3];
    Vector4f t;

    /// 3x3 sub-derminants 
    Vector4f v0( r2.x, r2.x, r3.x, r2.x );
    Vector4f v1( r2.y, r2.y, r3.y, r2.y );
    Vector4f v2( r2.z, r2.z, r3.z, r2.z );
    Vector4f v3( r2.w, r2.w, r3.w, r2.w );

    Vector4f w0( d25.x, d23.z, d22.x, d22.x );
    Vector4f w1( d25.y, d23.w, d22.y, d22.y );
    Vector4f w2( d24.z, d23.x, d21.z, d21.z );
    Vector4f w3( d24.w, d23.y, d21.w, d21.w );
    Vector4f w4( d24.y, d22.w, d21.y, d21.y );
    Vector4f w5( d24.x, d22.z, d21.x, d21.x );

    v[0]=v1; v[0]*=w1;   t=v2; t*=w0; v[0]-=t;   t=v3; t*=w3; v[0]+=t;
    v[1]=v0; v[1]*=w1;   t=v2; t*=w2; v[1]-=t;   t=v3; t*=w4; v[1]+=t;
    v[2]=v0; v[2]*=w0;   t=v1; t*=w2; v[2]-=t;   t=v3; t*=w5; v[2]+=t;
    v[3]=v0; v[3]*=w3;   t=v1; t*=w4; v[3]-=t;   t=v2; t*=w5; v[3]+=t;

    v[0].x = -v[0].x;
    v[0].w = -v[0].w;
    v[1].y = -v[1].y;
    v[1].z = -v[1].z;
    v[2].x = -v[2].x;
    v[2].w = -v[2].w;
    v[3].y = -v[3].y;
    v[3].z = -v[3].z;

    float inv_det = 1.0f / (
        m[ 3] * r3.x + 
        m[ 7] * r3.y + 
        m[11] * r3.z + 
        m[15] * r3.w
    );

    (*this) *= inv_det;
}

/// \brief Sets 'this' to the identity Matrix
void Matrix4f::identity(){
#if _SIMD_x86
    r0 = _mm_setr_ps(1.0, 0.0, 0.0, 0.0);
    r1 = _mm_setr_ps(0.0, 1.0, 0.0, 0.0);
    r2 = _mm_setr_ps(0.0, 0.0, 1.0, 0.0);
    r3 = _mm_setr_ps(0.0, 0.0, 0.0, 1.0);
#else 
    m[ 0]=1.0; m[ 1]=0.0; m[ 2]=0.0; m[ 3]=0.0;
    m[ 4]=0.0; m[ 5]=1.0; m[ 6]=0.0; m[ 7]=0.0;
    m[ 8]=0.0; m[ 9]=0.0; m[10]=1.0; m[11]=0.0;
    m[12]=0.0; m[13]=0.0; m[14]=0.0; m[15]=1.0;
#endif
}

/// \brief In-place transposition of 'this' Matrix
void Matrix4f::transpose(){
    float t;
    t=m[ 4]; m[ 4]=m[ 1]; m[ 1]=t;
    t=m[ 8]; m[ 8]=m[ 2]; m[ 2]=t;
    t=m[ 9]; m[ 9]=m[ 6]; m[ 6]=t;
    t=m[12]; m[12]=m[ 3]; m[ 3]=t;
    t=m[13]; m[13]=m[ 7]; m[ 7]=t;
    t=m[14]; m[14]=m[11]; m[11]=t;
}

/////////////////////////////////////////////////////////////////
// double 
/////////////////////////////////////////////////////////////////
Matrix4d::Matrix4d(){
#ifdef _SIMD_x86
    r0 = _mm256_setzero_pd();
    r1 = _mm256_setzero_pd();
    r2 = _mm256_setzero_pd();
    r3 = _mm256_setzero_pd();
#else 
    m[ 0]=0.0; m[ 1]=0.0; m[ 2]=0.0; m[ 3]=0.0;
    m[ 4]=0.0; m[ 5]=0.0; m[ 6]=0.0; m[ 7]=0.0;
    m[ 8]=0.0; m[ 9]=0.0; m[10]=0.0; m[11]=0.0;
    m[12]=0.0; m[13]=0.0; m[14]=0.0; m[15]=0.0;
#endif 
}

Matrix4d::Matrix4d(const Matrix4d &A){
#ifdef _SIMD_x86
    r0 = A.r0;
    r1 = A.r1;
    r2 = A.r2;
    r3 = A.r3;
#else 
    m[ 0]=A.m[ 0]; m[ 1]=A.m[ 1]; m[ 2]=A.m[ 2]; m[ 3]=A.m[ 3];
    m[ 4]=A.m[ 4]; m[ 5]=A.m[ 5]; m[ 6]=A.m[ 6]; m[ 7]=A.m[ 7];
    m[ 8]=A.m[ 8]; m[ 9]=A.m[ 9]; m[10]=A.m[10]; m[11]=A.m[11];
    m[12]=A.m[12]; m[13]=A.m[13]; m[14]=A.m[14]; m[15]=A.m[15];
#endif 
}

Matrix4d &Matrix4d::operator = (const Matrix4d &A){
#ifdef _SIMD_x86
    r0 = A.r0;
    r1 = A.r1;
    r2 = A.r2;
    r3 = A.r3;
#else 
    m[ 0]=A.m[ 0]; m[ 1]=A.m[ 1]; m[ 2]=A.m[ 2]; m[ 3]=A.m[ 3];
    m[ 4]=A.m[ 4]; m[ 5]=A.m[ 5]; m[ 6]=A.m[ 6]; m[ 7]=A.m[ 7];
    m[ 8]=A.m[ 8]; m[ 9]=A.m[ 9]; m[10]=A.m[10]; m[11]=A.m[11];
    m[12]=A.m[12]; m[13]=A.m[13]; m[14]=A.m[14]; m[15]=A.m[15];
#endif 
    return *this;
}

Matrix4d& Matrix4d::operator = (double s){
#ifdef _SIMD_x86
    r0 = _mm256_set1_pd(s);
    r1 = _mm256_set1_pd(s);
    r2 = _mm256_set1_pd(s);
    r3 = _mm256_set1_pd(s);
#else 
    m[ 0]=s; m[ 1]=s; m[ 2]=s; m[ 3]=s;
    m[ 4]=s; m[ 5]=s; m[ 6]=s; m[ 7]=s;
    m[ 8]=s; m[ 9]=s; m[10]=s; m[11]=s;
    m[12]=s; m[13]=s; m[14]=s; m[15]=s;
#endif 
    return *this;
}

Matrix4d& Matrix4d::operator += (const Matrix4d &A){
#ifdef _SIMD_x86
    r0 = _mm256_add_pd(r0, A.r0);
    r1 = _mm256_add_pd(r1, A.r1);
    r2 = _mm256_add_pd(r2, A.r2);
    r3 = _mm256_add_pd(r3, A.r3);
#else 
    m[ 0]+=A.m[ 0]; m[ 1]+=A.m[ 1]; m[ 2]+=A.m[ 2]; m[ 3]+=A.m[ 3];
    m[ 4]+=A.m[ 4]; m[ 5]+=A.m[ 5]; m[ 6]+=A.m[ 6]; m[ 7]+=A.m[ 7];
    m[ 8]+=A.m[ 8]; m[ 9]+=A.m[ 9]; m[10]+=A.m[10]; m[11]+=A.m[11];
    m[12]+=A.m[12]; m[13]+=A.m[13]; m[14]+=A.m[14]; m[15]+=A.m[15];
#endif 
    return *this;
}

Matrix4d& Matrix4d::operator -= (const Matrix4d &A){
#ifdef _SIMD_x86
    r0 = _mm256_sub_pd(r0, A.r0);
    r1 = _mm256_sub_pd(r1, A.r1);
    r2 = _mm256_sub_pd(r2, A.r2);
    r3 = _mm256_sub_pd(r3, A.r3);
#else 
    m[ 0]-=A.m[ 0]; m[ 1]-=A.m[ 1]; m[ 2]-=A.m[ 2]; m[ 3]-=A.m[ 3];
    m[ 4]-=A.m[ 4]; m[ 5]-=A.m[ 5]; m[ 6]-=A.m[ 6]; m[ 7]-=A.m[ 7];
    m[ 8]-=A.m[ 8]; m[ 9]-=A.m[ 9]; m[10]-=A.m[10]; m[11]-=A.m[11];
    m[12]-=A.m[12]; m[13]-=A.m[13]; m[14]-=A.m[14]; m[15]-=A.m[15];
#endif 
    return *this;
}

Matrix4d& Matrix4d::operator *= (const Matrix4d &A){
#ifdef _SIMD_x86
    __m256d row;
    Vector4d t;

    t.r = r0;
    row = _mm256_mul_pd  (A.r0, _mm256_set1_pd(t.x));
    row = _mm256_fmadd_pd(A.r1, _mm256_set1_pd(t.y), row);
    row = _mm256_fmadd_pd(A.r2, _mm256_set1_pd(t.z), row);
    row = _mm256_fmadd_pd(A.r3, _mm256_set1_pd(t.w), row);
    r0 = row;

    t.r = r1;
    row = _mm256_mul_pd  (A.r0, _mm256_set1_pd(t.x));
    row = _mm256_fmadd_pd(A.r1, _mm256_set1_pd(t.y), row);
    row = _mm256_fmadd_pd(A.r2, _mm256_set1_pd(t.z), row);
    row = _mm256_fmadd_pd(A.r3, _mm256_set1_pd(t.w), row);
    r1 = row;

    t.r = r2;
    row = _mm256_mul_pd  (A.r0, _mm256_set1_pd(t.x));
    row = _mm256_fmadd_pd(A.r1, _mm256_set1_pd(t.y), row);
    row = _mm256_fmadd_pd(A.r2, _mm256_set1_pd(t.z), row);
    row = _mm256_fmadd_pd(A.r3, _mm256_set1_pd(t.w), row);
    r2 = row;
    
    t.r = r3;
    row = _mm256_mul_pd  (A.r0, _mm256_set1_pd(t.x));
    row = _mm256_fmadd_pd(A.r1, _mm256_set1_pd(t.y), row);
    row = _mm256_fmadd_pd(A.r2, _mm256_set1_pd(t.z), row);
    row = _mm256_fmadd_pd(A.r3, _mm256_set1_pd(t.w), row);
    r3 = row;
#else 
    double t[4];

    t[0]=m[ 0]; t[1]=m[ 1]; t[2]=m[ 2]; t[3]=m[ 3];
    m[ 0]=t[0]*A.m[ 0] + t[1]*A.m[ 4] + t[2]*A.m[ 8] + t[3]*A.m[12];
    m[ 1]=t[0]*A.m[ 1] + t[1]*A.m[ 5] + t[2]*A.m[ 9] + t[3]*A.m[13];
    m[ 2]=t[0]*A.m[ 2] + t[1]*A.m[ 6] + t[2]*A.m[10] + t[3]*A.m[14];
    m[ 3]=t[0]*A.m[ 3] + t[1]*A.m[ 7] + t[2]*A.m[11] + t[3]*A.m[15];

    t[0]=m[ 4]; t[1]=m[ 5]; t[2]=m[ 6]; t[3]=m[ 7];
    m[ 4]=t[0]*A.m[ 0] + t[1]*A.m[ 4] + t[2]*A.m[ 8] + t[3]*A.m[12];
    m[ 5]=t[0]*A.m[ 1] + t[1]*A.m[ 5] + t[2]*A.m[ 9] + t[3]*A.m[13];
    m[ 6]=t[0]*A.m[ 2] + t[1]*A.m[ 6] + t[2]*A.m[10] + t[3]*A.m[14];
    m[ 7]=t[0]*A.m[ 3] + t[1]*A.m[ 7] + t[2]*A.m[11] + t[3]*A.m[15];

    t[0]=m[ 8]; t[1]=m[ 9]; t[2]=m[10]; t[3]=m[11];
    m[ 8]=t[0]*A.m[ 0] + t[1]*A.m[ 4] + t[2]*A.m[ 8] + t[3]*A.m[12];
    m[ 9]=t[0]*A.m[ 1] + t[1]*A.m[ 5] + t[2]*A.m[ 9] + t[3]*A.m[13];
    m[10]=t[0]*A.m[ 2] + t[1]*A.m[ 6] + t[2]*A.m[10] + t[3]*A.m[14];
    m[11]=t[0]*A.m[ 3] + t[1]*A.m[ 7] + t[2]*A.m[11] + t[3]*A.m[15];

    t[0]=m[12]; t[1]=m[13]; t[2]=m[14]; t[3]=m[15];
    m[12]=t[0]*A.m[ 0] + t[1]*A.m[ 4] + t[2]*A.m[ 8] + t[3]*A.m[12];
    m[13]=t[0]*A.m[ 1] + t[1]*A.m[ 5] + t[2]*A.m[ 9] + t[3]*A.m[13];
    m[14]=t[0]*A.m[ 2] + t[1]*A.m[ 6] + t[2]*A.m[10] + t[3]*A.m[14];
    m[15]=t[0]*A.m[ 3] + t[1]*A.m[ 7] + t[2]*A.m[11] + t[3]*A.m[15];
#endif
    return *this;
}

Matrix4d& Matrix4d::operator *= (double s){
#ifdef _SIMD_x86
    __m256d vs = _mm256_set1_pd(s);
    r0 = _mm256_mul_pd(r0, vs);
    r1 = _mm256_mul_pd(r1, vs);
    r2 = _mm256_mul_pd(r2, vs);
    r3 = _mm256_mul_pd(r3, vs);
#else 
    m[ 0]*=s; m[ 1]*=s; m[ 2]*=s; m[ 3]*=s;
    m[ 4]*=s; m[ 5]*=s; m[ 6]*=s; m[ 7]*=s;
    m[ 8]*=s; m[ 9]*=s; m[10]*=s; m[11]*=s;
    m[12]*=s; m[13]*=s; m[14]*=s; m[15]*=s;
#endif 
    return *this;
}

/// \brief Sets 'this' to the identity Matrix
void Matrix4d::identity(){
#if _SIMD_x86
    r0 = _mm256_setr_pd(1.0, 0.0, 0.0, 0.0);
    r1 = _mm256_setr_pd(0.0, 1.0, 0.0, 0.0);
    r2 = _mm256_setr_pd(0.0, 0.0, 1.0, 0.0);
    r3 = _mm256_setr_pd(0.0, 0.0, 0.0, 1.0);
#else 
    m[ 0]=1.0; m[ 1]=0.0; m[ 2]=0.0; m[ 3]=0.0;
    m[ 4]=0.0; m[ 5]=1.0; m[ 6]=0.0; m[ 7]=0.0;
    m[ 8]=0.0; m[ 9]=0.0; m[10]=1.0; m[11]=0.0;
    m[12]=0.0; m[13]=0.0; m[14]=0.0; m[15]=1.0;
#endif
}

/// \brief In-place transposition of 'this' Matrix
void Matrix4d::transpose(){
    double t;
    t=m[ 4]; m[ 4]=m[ 1]; m[ 1]=t;
    t=m[ 8]; m[ 8]=m[ 2]; m[ 2]=t;
    t=m[ 9]; m[ 9]=m[ 6]; m[ 6]=t;
    t=m[12]; m[12]=m[ 3]; m[ 3]=t;
    t=m[13]; m[13]=m[ 7]; m[ 7]=t;
    t=m[14]; m[14]=m[11]; m[11]=t;
}

/// \brief In-place computation of the 4x4 inverse Matrix
/// \note \f$ A*A^-1 = A^-1*A = Identity \f$
/// \attention Doesn't check if Matrix is singular
void Matrix4d::invert()
{
    Vector4d tu;
    Vector4d tv;

    /// 2x2 sub-derminants
    Vector4d d21;
    tu = Vector4d(m[ 0], m[ 0], m[ 0], m[ 1]);
    tv = Vector4d(m[ 5], m[ 6], m[ 7], m[ 6]);
    tu *= tv;
    d21 = tu;
    tu = Vector4d(m[ 1], m[ 2], m[ 3], m[ 2]);
    tv = Vector4d(m[ 4], m[ 4], m[ 4], m[ 5]);
    tu *= tv;
    d21 -= tu;
    
    Vector4d d22;
    tu = Vector4d(m[ 1], m[ 2], m[ 0], m[ 0]);
    tv = Vector4d(m[ 7], m[ 7], m[13], m[14]);
    tu *= tv;
    d22 = tu;
    tu = Vector4d(m[ 3], m[ 3], m[ 1], m[ 2]);
    tv = Vector4d(m[ 5], m[ 6], m[12], m[12]);
    tu *= tv;
    d22 -= tu;
    
    Vector4d d23;
    tu = Vector4d(m[ 0], m[ 1], m[ 1], m[ 2]);
    tv = Vector4d(m[15], m[14], m[15], m[15]);
    tu *= tv;
    d23 = tu;
    tu = Vector4d(m[ 3], m[ 2], m[ 3], m[ 3]);
    tv = Vector4d(m[12], m[13], m[13], m[14]);
    tu *= tv;
    d23 -= tu;
    
    Vector4d d24;
    tu = Vector4d(m[ 4], m[ 4], m[ 4], m[ 5]);
    tv = Vector4d(m[13], m[14], m[15], m[14]);
    tu *= tv;
    d24 = tu;
    tu = Vector4d(m[ 5], m[ 6], m[ 7], m[ 6]);
    tv = Vector4d(m[12], m[12], m[12], m[13]);
    tu *= tv;
    d24 -= tu;
    
    Vector4d d25;
    tu = Vector4d(m[ 5], m[ 6], 0, 0);
    tv = Vector4d(m[15], m[15], 0, 0);
    tu *= tv;
    d25 = tu;
    tu = Vector4d(m[ 7], m[ 7], 0, 0);
    tv = Vector4d(m[13], m[14], 0, 0);
    tu *= tv;
    d25 -= tu;

    /// temp 
    Vector4d r2 = v[2];
    Vector4d r3 = v[3];
    Vector4d t;

    /// 3x3 sub-derminants 
    Vector4d v0( r2.x, r2.x, r3.x, r2.x );
    Vector4d v1( r2.y, r2.y, r3.y, r2.y );
    Vector4d v2( r2.z, r2.z, r3.z, r2.z );
    Vector4d v3( r2.w, r2.w, r3.w, r2.w );

    Vector4d w0( d25.x, d23.z, d22.x, d22.x );
    Vector4d w1( d25.y, d23.w, d22.y, d22.y );
    Vector4d w2( d24.z, d23.x, d21.z, d21.z );
    Vector4d w3( d24.w, d23.y, d21.w, d21.w );
    Vector4d w4( d24.y, d22.w, d21.y, d21.y );
    Vector4d w5( d24.x, d22.z, d21.x, d21.x );

    v[0]=v1; v[0]*=w1;   t=v2; t*=w0; v[0]-=t;   t=v3; t*=w3; v[0]+=t;
    v[1]=v0; v[1]*=w1;   t=v2; t*=w2; v[1]-=t;   t=v3; t*=w4; v[1]+=t;
    v[2]=v0; v[2]*=w0;   t=v1; t*=w2; v[2]-=t;   t=v3; t*=w5; v[2]+=t;
    v[3]=v0; v[3]*=w3;   t=v1; t*=w4; v[3]-=t;   t=v2; t*=w5; v[3]+=t;

    v[0].x = -v[0].x;
    v[0].w = -v[0].w;
    v[1].y = -v[1].y;
    v[1].z = -v[1].z;
    v[2].x = -v[2].x;
    v[2].w = -v[2].w;
    v[3].y = -v[3].y;
    v[3].z = -v[3].z;

    double inv_det = 1.0 / (
        m[ 3] * r3.x + 
        m[ 7] * r3.y + 
        m[11] * r3.z + 
        m[15] * r3.w
    );

    (*this) *= inv_det;
}

/// \brief Uses this Matrix to transform the Vector \p b such that: \f$ a=M*b \f$
/// \param a (dst) Transformed Vector
/// \param b (src) Vector to be transformed
void Matrix4d::transform(Vector4d &a, const Vector4d &b){
#ifdef _SIMD_x86
    __m256d k;
    __m256d col0 = _mm256_setr_pd(m[0], m[4], m[ 8], m[12]);
    __m256d col1 = _mm256_setr_pd(m[1], m[5], m[ 9], m[13]);
    __m256d col2 = _mm256_setr_pd(m[2], m[6], m[10], m[14]);
    __m256d col3 = _mm256_setr_pd(m[3], m[7], m[11], m[15]);

    k = _mm256_mul_pd  (col0, _mm256_set1_pd(b.x));
    k = _mm256_fmadd_pd(col1, _mm256_set1_pd(b.y), k);
    k = _mm256_fmadd_pd(col2, _mm256_set1_pd(b.z), k);
    k = _mm256_fmadd_pd(col3, _mm256_set1_pd(b.w), k);

    a.r = k;
#else
    a.v[0] = m[ 0]*b.v[0] + m[ 1]*b.v[1] + m[ 2]*b.v[2] + m[ 3]*b.v[3];
    a.v[1] = m[ 4]*b.v[0] + m[ 5]*b.v[1] + m[ 6]*b.v[2] + m[ 7]*b.v[3];
    a.v[2] = m[ 8]*b.v[0] + m[ 9]*b.v[1] + m[10]*b.v[2] + m[11]*b.v[3];
    a.v[3] = m[12]*b.v[0] + m[13]*b.v[1] + m[14]*b.v[2] + m[15]*b.v[3];
#endif
}

/// Calculates the eigenvalues of a symmetric 3x3 matrix using Cardano's
/// analytical algorithm.
void Matrix4d::eigval3(Vector4d &ev) const {
    /// Determine coefficients of characteristic poynomial. 
    ///       | a   d   f |
    ///  A =  | d*  b   e |
    ///       | f*  e*  c |
    double ab = m[0]*m[5];
    double de = m[1]*m[6];
    double dd = m[1]*m[1];
    double ee = m[6]*m[6];
    double ff = m[2]*m[2];

    double tr = m[0] +m[5] +m[10];
    
    /// a*b + a*c + b*c - d^2 - e^2 - f^2
    double c1 = ab +m[0]*m[10] +m[5]*m[10] -dd -ee -ff; 
    
    /// c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e
    double c0 = m[10]*dd +m[0]*ee +m[5]*ff -m[10]*ab -m[2]*de*2.0;

    double p = tr*tr -3.0*c1;
    double q = tr*(p -1.5*c1) -13.5*c0;
    double r = std::sqrt( std::fabs(p));

    double phi;
    phi = 27.0* (0.25*(c1*c1) * (p-c1) + c0*(q +6.75*c0));
    phi = std::atan2(std::sqrt(std::fabs(phi)), q)*(1.0/3.0);

    double c = r * std::cos(phi);
    double s = r * std::sin(phi) * (1.0/std::sqrt(3.0));

    ev.y = (tr-c)*(1.0/3.0);
    ev.z = ev.y + s;
    ev.x = ev.y + c;
    ev.y -= s;

    /// \todo Sort by absolute value? change below
    ///if (ev.x < ev.z) std::swap(ev.x, ev.z);
    ///if (ev.x < ev.y) std::swap(ev.x, ev.y);
    ///if (ev.y < ev.z) std::swap(ev.y, ev.z);
}

std::ostream &operator << (std::ostream &os, const Vector4f &a){
    os << std::fixed;
    os << std::setprecision(6);
    os << std::setw(12) << a.x << ' ';
    os << std::setw(12) << a.y << ' ';
    os << std::setw(12) << a.z << ' ';
    os << std::setw(12) << a.w;
    return os;
}
std::ostream &operator << (std::ostream &os, const Vector4d &a){
    os << std::fixed;
    os << std::setprecision(6);
    os << std::setw(12) << a.x << ' ';
    os << std::setw(12) << a.y << ' ';
    os << std::setw(12) << a.z << ' ';
    os << std::setw(12) << a.w;
    return os;
}
std::ostream &operator << (std::ostream &os, const Matrix4f &a) {
    os << std::fixed;
    os << std::setprecision(6);
    for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
        os << std::setw(12);
        os << a.m[i*4 + j] << ' ';
        }

        if (i < 3) 
            os << std::endl;
    }
    return os;
}
std::ostream &operator << (std::ostream &os, const Matrix4d &a) {
    os << std::fixed;
    os << std::setprecision(6);
    for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
        os << std::setw(12);
        os << a.m[i*4 + j] << ' ';
        }

        if (i < 3) 
            os << std::endl;
    }
    return os;
}

}/// !namespace math
