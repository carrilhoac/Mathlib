
#ifndef _Math_Functions_h
#define _Math_Functions_h 1

////////////////////////////////////////////////////////////////////////////
// 3D math
template<typename TA, typename TB>
inline MatrixGeneric<TA> 
cross(const MatrixGeneric<TA>& a, const MatrixGeneric<TB>& b)
{
  MatrixGeneric<TA> c(3);
  c(0) = a(1) * b(2) - a(2) * b(1);
  c(1) = a(2) * b(0) - a(0) * b(2);
  c(2) = a(0) * b(1) - a(1) * b(0);
  return c;
}

template<typename TA, typename TB, typename TC>
inline MatrixGeneric<TA> 
normal(
  const MatrixGeneric<TA>& v1, 
  const MatrixGeneric<TB>& v2, 
  const MatrixGeneric<TC>& v3)
{
  return unit(cross(v2-v1, v3-v1));
}

template<typename TA>
inline MatrixGeneric<TA> 
rotationX(const typename TA::DataType &thetaRad) 
{
  using Scalar = typename TA::DataType;
  const Scalar thetaSin = std::sin(thetaRad);
  const Scalar thetaCos = std::cos(thetaRad);

  MatrixGeneric<TA> R(3,3);
  R.setIdentity();

  R(1,1) =  thetaCos;
  R(2,2) =  thetaCos;
  R(1,2) = -thetaSin;
  R(2,1) =  thetaSin;
  return R;
}

template<typename TA>
inline MatrixGeneric<TA> 
rotationY(const typename TA::DataType &thetaRad) 
{
  using Scalar = typename TA::DataType;
  const Scalar thetaSin = std::sin(thetaRad);
  const Scalar thetaCos = std::cos(thetaRad);

  MatrixGeneric<TA> R(3,3);
  R.setIdentity();

  R(0,0) =  thetaCos;
  R(2,2) =  thetaCos;
  R(0,2) =  thetaSin;
  R(2,0) = -thetaSin;
  return R;
}

template<typename TA>
inline MatrixGeneric<TA> 
rotationZ(const typename TA::DataType &thetaRad) 
{
  using Scalar = typename TA::DataType;
  const Scalar thetaSin = std::sin(thetaRad);
  const Scalar thetaCos = std::cos(thetaRad);

  MatrixGeneric<TA> R(3,3);
  R.setIdentity();

  R(0,0) =  thetaCos;
  R(1,1) =  thetaCos;
  R(0,1) = -thetaSin;
  R(1,0) =  thetaSin;
  return R;
}
////////////////////////////////////////////////////////////////////////////

template<typename MatrixLayout>
inline MatrixGeneric<MatrixLayout> 
transpose(const MatrixGeneric<MatrixLayout>& A) 
{
  const int m = A.rows();
  const int n = A.cols();

  MatrixGeneric<MatrixLayout> T(n, m);
  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    T(j,i) = A(i,j);

  return T;
}

template<typename MatrixLayout>
inline MatrixGeneric<MatrixLayout> 
inv(const MatrixGeneric<MatrixLayout>& A) 
{
  return LU(A).inverse();
}

template<typename MatrixLayout>
inline typename MatrixLayout::DataType 
det(const MatrixGeneric<MatrixLayout>& A) 
{
  return LU(A).det();
}

template<typename MatrixLayout>
inline MatrixGeneric<MatrixLayout> 
nullspace(const MatrixGeneric<MatrixLayout>& A) 
{
  return SVD(A).nullspace();
}

template<typename MatrixLayout>
inline int 
rank(const MatrixGeneric<MatrixLayout>& A, 
  typename MatrixLayout::DataType tol = eps<typename MatrixLayout::DataType>()) 
{
  return SVD(A).rank(tol);
}

template<typename MatrixLayout>
inline typename MatrixLayout::DataType 
cond(const MatrixGeneric<MatrixLayout>& A) 
{
  return SVD(A).cond(); // same equation as in MATLAB, scilab, octave...
}

template<typename MatrixLayout>
inline typename MatrixLayout::DataType 
norm2(const MatrixGeneric<MatrixLayout>& A) 
{
  return SVD(A).norm2();
}

template<typename MatrixLayout>
inline typename MatrixLayout::DataType 
trace(const MatrixGeneric<MatrixLayout>& A) 
{
  using Scalar = typename MatrixLayout::DataType;
  const int n = A.rows();
  Scalar trace_(0);

  for (int i = 0; i < n; i++)
    trace_ += A(i,i);

  return trace_;
}

template<typename TA, typename TB>
inline typename TA::DataType 
dot(const MatrixGeneric<TA>& a, const MatrixGeneric<TB>& b)
{
  using Scalar = typename TA::DataType;
  const int n = a.rows();
  Scalar dot_(0);

  for (int i = 0; i < n; i++)
    dot_ += a(i) * b(i);

  return dot_;
}

template<typename TA>
inline typename TA::DataType 
norm(const MatrixGeneric<TA>& a)
{
  using Scalar = typename TA::DataType;
  const int m = a.rows();
  const int n = a.cols();
  Scalar norm_(0);

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    norm_ += a(i,j) * a(i,j);

  return std::sqrt(norm_);
}

template<typename TA>
inline MatrixGeneric<TA> 
unit(const MatrixGeneric<TA>& a)
{
  using Scalar = typename TA::DataType;
  return a * (Scalar(1) / norm(a)); 
}

template<typename TA, typename TB>
inline typename TA::DataType 
distanceSqr(const MatrixGeneric<TA>& a, const MatrixGeneric<TB>& b)
{
  using Scalar = typename TA::DataType;
  const int n = a.rows();
  Scalar distSqr_(0);

  for (int i = 0; i < n; i++) {
    const Scalar delta_ = a(i) - b(i);
    distSqr_ += delta_ * delta_;
  }

  return distSqr_;
}

template<typename TA, typename TB>
inline typename TA::DataType 
distance(const MatrixGeneric<TA>& a, const MatrixGeneric<TB>& b)
{
  return std::sqrt(distanceSqr(a,b));
}

#endif 

