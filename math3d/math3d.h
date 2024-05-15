
/////////////////////////////////////////////////////////////////
/// \author Andre Caceres Carrilho
/// \namespace math  
/// \file math3d.h
/// \todo implement NEON SIMD for ARM (float32 only)
/////////////////////////////////////////////////////////////////

#ifndef _Math3D_h
#define _Math3D_h 1

#if (__cplusplus < 201703L)
    #error C++17 or later is required 
    /// C++17 or later is needed in order to have proper alignment of 
    /// SIMD data inside std::vector<> containers (otherwise segfault)
#endif
#if defined(__AVX2__)
    #include <immintrin.h>
    #define _SIMD_x86 1 /// Assuming we have both AVX2 and FMA3 
#endif 
#if defined(__ARM_NEON__)
    #include <arm_neon.h>
    #define _SIMD_ARM 1
#endif 

#include <iostream>
#include <iomanip>

/// \brief 3D math classes and methods
/// \note Supports AVX2 and FMA on x86_64 CPUs. ARM NEON only for float.
namespace math{

/////////////////////////////////////////////////////////////////
/// \class Vector4f
/// \brief \f$ Vector4f (x,y,z,w) \f$. 
/// \note Data type: float (f32).
/////////////////////////////////////////////////////////////////
class alignas(16) Vector4f {
public:
    union {
        struct { float x, y, z, w; };
        float v[4];
#ifdef _SIMD_x86
        __m128 r;
#endif 
#ifdef _SIMD_ARM
        float32x4_t r;
#endif 
    };
public:
    Vector4f(float X=0.0f, float Y=0.0f, float Z=0.0f, float W=0.0f);
    Vector4f(const Vector4f &a);
    Vector4f& operator = (const Vector4f &a);  
    Vector4f& operator = (float s);
    Vector4f& operator += (const Vector4f &a); 
    Vector4f& operator -= (const Vector4f &a);
    Vector4f& operator *= (const Vector4f &a);
    Vector4f& operator *= (float s);
    float norm() const;
    void normalize();
    void homogenize();
    Vector4f& max(const Vector4f &a);
    Vector4f& min(const Vector4f &a);
    bool operator > (const Vector4f &a) const;
    bool operator < (const Vector4f &a) const;
    bool operator >= (const Vector4f &a) const;
    bool operator <= (const Vector4f &a) const;
    friend std::ostream &operator << (std::ostream &os, const Vector4f &a);
};

/////////////////////////////////////////////////////////////////
/// \class Matrix4f 
/// \brief 4x4 Matrix (row-major).
/// \note Data type: float (f32).
/////////////////////////////////////////////////////////////////
class alignas(64) Matrix4f {
public:
    union {
        float m[16];
        Vector4f v[4];
#ifdef _SIMD_x86
        struct { __m128 r0, r1, r2, r3; };
#endif 
#ifdef _SIMD_ARM
        struct { float32x4_t r0, r1, r2, r3; };
#endif 
    };
public:
    Matrix4f();
    Matrix4f(const Matrix4f &A);
    Matrix4f& operator = (const Matrix4f &A);
    Matrix4f& operator = (float s);
    Matrix4f& operator += (const Matrix4f &A);
    Matrix4f& operator -= (const Matrix4f &A);
    Matrix4f& operator *= (const Matrix4f &A);
    Matrix4f& operator *= (float s);
    void identity();
    void transpose();
    void invert();
    void transform(Vector4f &a, const Vector4f &b);
    friend std::ostream &operator << (std::ostream &os, const Matrix4f &a);
};

/////////////////////////////////////////////////////////////////
/// \class Vector4d 
/// \brief \f$ Vector4d (x,y,z,w) \f$. 
/// \note Data type: double (f64).
/////////////////////////////////////////////////////////////////
class alignas(32) Vector4d {
public:
    union {
        struct { double x, y, z, w; };
        double v[4];
#ifdef _SIMD_x86
        __m256d r;
#endif 
    };
public:
    Vector4d(double X=0.0, double Y=0.0, double Z=0.0, double W=0.0);
    Vector4d(const Vector4d &a);
    Vector4d& operator = (const Vector4d &a);  
    Vector4d& operator = (double s);      
    Vector4d& operator += (const Vector4d &a); 
    Vector4d& operator -= (const Vector4d &a);
    Vector4d& operator *= (const Vector4d &a);
    Vector4d& operator *= (double s);
    double norm() const;
    void normalize();
    void homogenize();
    Vector4d& max(const Vector4d &a);
    Vector4d& min(const Vector4d &a);
    bool operator > (const Vector4d &a) const;
    bool operator < (const Vector4d &a) const;
    bool operator >= (const Vector4d &a) const;
    bool operator <= (const Vector4d &a) const;
    friend std::ostream &operator << (std::ostream &os, const Vector4d &a);
};

/////////////////////////////////////////////////////////////////
/// \class Matrix4d 
/// \brief 4x4 Matrix (row-major).
/// \note Data type: double (f64).
/////////////////////////////////////////////////////////////////
class alignas(128) Matrix4d {
public:
    union {
        double m[16];
        Vector4d v[4];
#ifdef _SIMD_x86
        struct { __m256d r0, r1, r2, r3; };
#endif 
    };
public:
    Matrix4d();
    Matrix4d(const Matrix4d &A);
    Matrix4d& operator = (const Matrix4d &A);
    Matrix4d& operator = (double s);
    Matrix4d& operator += (const Matrix4d &A);
    Matrix4d& operator -= (const Matrix4d &A);
    Matrix4d& operator *= (const Matrix4d &A);
    Matrix4d& operator *= (double s);
    void identity();
    void transpose();
    void invert();
    void transform(Vector4d &a, const Vector4d &b);
    void eigval3(Vector4d &ev) const;
    friend std::ostream &operator << (std::ostream &os, const Matrix4d &a);
};

double dist2(const Vector4d &a, const Vector4d &b);
float dist2(const Vector4f &a, const Vector4f &b);
double dist(const Vector4d &a, const Vector4d &b);
float dist(const Vector4f &a, const Vector4f &b);
double dot(const Vector4d &a, const Vector4d &b);
float dot(const Vector4f &a, const Vector4f &b);
void cross(Vector4d &c, const Vector4d &a, const Vector4d &b);
void cross(Vector4f &c, const Vector4f &a, const Vector4f &b);
void corr(Matrix4d &C, const Vector4d &a, const Vector4d &b);
void corr(Matrix4f &C, const Vector4f &a, const Vector4f &b);
void corr(Matrix4d &C, const Vector4d &a);  /// https://stackoverflow.com/a/38324464
void corr(Matrix4f &C, const Vector4f &a);

}/// !namespace math
#endif /// !_Math3D_h
