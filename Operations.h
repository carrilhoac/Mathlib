
#ifndef _Math_Operations_h
#define _Math_Operations_h 1

template<typename TA, typename TB>
inline MatrixGeneric<TA>& operator += (
  MatrixGeneric<TA>& A, const MatrixGeneric<TB>& B) 
{
  int m = A.rows();
  int n = A.cols();

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    A(i,j) += B(i,j);

  return A;
}

template<typename TA, typename TB>
inline MatrixGeneric<TA>& operator -= (
  MatrixGeneric<TA>& A, const MatrixGeneric<TB>& B) 
{
  int m = A.rows();
  int n = A.cols();

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    A(i,j) -= B(i,j);

  return A;
}

template<typename TA, typename TB>
inline MatrixGeneric<TA> operator + (
  const MatrixGeneric<TA>& A, const MatrixGeneric<TB>& B) 
{
  MatrixGeneric<TA> C(A);
  C += B;
  return C;
}

template<typename TA, typename TB>
inline MatrixGeneric<TA> operator - (
  const MatrixGeneric<TA>& A, const MatrixGeneric<TB>& B) 
{
  MatrixGeneric<TA> C(A);
  C -= B;
  return C;
}

template<typename TA, typename TB>
inline MatrixGeneric<TA> operator * (
  const MatrixGeneric<TA>& A, const MatrixGeneric<TB>& B) 
{
  int m = A.rows();
  int n = A.cols(); // must be equal to B.rows()
  int q = B.cols();

  MatrixGeneric<TA> C(m, q);
  C.setZero();

  for (int i = 0; i < m; i++)
  for (int k = 0; k < n; k++)
  for (int j = 0; j < q; j++)
    C(i,j) += A(i,k) * B(k,j);

  return C;
}

template<typename TA, typename Scalar>
inline MatrixGeneric<TA> operator * (
  const MatrixGeneric<TA>& A, Scalar s)
{
  typename MatrixGeneric<TA>::DataType x(s);
  MatrixGeneric<TA> B(A);

  int m = B.rows();
  int n = B.cols();

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    B(i,j) *= x;

  return B;
}

template<typename TA, typename Scalar>
inline MatrixGeneric<TA>& operator *= (
  MatrixGeneric<TA>& A, Scalar s)
{
  typename MatrixGeneric<TA>::DataType x(s);

  int m = A.rows();
  int n = A.cols();

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    A(i,j) *= x;

  return A;
}

template<typename TA, typename Scalar>
inline MatrixGeneric<TA> operator / (
  const MatrixGeneric<TA>& A, Scalar s)
{
  typename MatrixGeneric<TA>::DataType x(1);
  return A * (x / s);
}

template<typename TA, typename Scalar>
inline MatrixGeneric<TA>& operator /= (
  MatrixGeneric<TA>& A, Scalar s)
{
  typename MatrixGeneric<TA>::DataType x(1);
  x /= s;
  A *= x;
  return A;
}

#endif
