
#ifndef _Math_Cholesky_h
#define _Math_Cholesky_h 1

// Cholesky factorization: A = L * L^T where  L is a lower triangular matrix.
// Only the lower triangular of A is accessed during the factorization.
// Matrix A must be symmetric semi-positive definite.
template<typename MatrixType>
class Cholesky {
protected:
  typedef typename MatrixType::DataType DataType;
  typedef typename std::vector<DataType> VectorType;

protected:
  MatrixType m_matrix;
  bool m_isPositiveDefinite;
  
public:
  Cholesky(const MatrixType& A)
    : m_matrix(A)
    , m_isPositiveDefinite(false)
  { factorize(); }

  bool isPositiveDefinite() const { return m_isPositiveDefinite; }
  MatrixType getL() const { return m_matrix; }
  MatrixType inverse() const;
  MatrixType solve(const MatrixType& b) const;

private:
  void factorize();
};

template<typename MatrixType>
inline void Cholesky<MatrixType>::factorize()
{
  if (! m_matrix.isSquare()) 
    return;

  int n = m_matrix.rows();
  VectorType m_invSqrt(n);
  DataType x(0);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      x = m_matrix(i,j);

      for (int k = 0; k < j; k++) 
        x -= m_matrix(i,k) * m_matrix(j,k);
    
      m_matrix(i,j) = x * m_invSqrt[j];
    }

    x = m_matrix(i,i);

    for (int k = 0; k < i; k++)
      x -= m_matrix(i,k) * m_matrix(i,k);

    if (x <= DataType(0))
      return;

    m_invSqrt[i] = DataType(1) / std::sqrt(x);
    m_matrix(i,i) = x * m_invSqrt[i];
  }

  m_isPositiveDefinite = true; 
  m_matrix.clearUpperTri();
}


// Special backsubstitution used for inverse computation
// Solve A*x = b , x = b being a row in the inverse matrix.
// The inverse is also symmetric, since inv(A) = transpose(inv(A)) 
// if A is S.P.D. so we only need to solve the triangular 

template<typename MatrixType>
inline MatrixType Cholesky<MatrixType>::inverse() const
{
  int n = m_matrix.rows();
  MatrixType B(n, n);

  // I think this should be checked before calling this method
  //if (! m_isPositiveDefinite)
  //  return B;

  B.setIdentity();

  for (int j = 0; j < n; j++) {

    //  Solve L*y = b
    for (int k = j; k < n; k++) {
      for (int i = j; i < k; i++) 
        B(j,k) -= B(j,i) * m_matrix(k,i);

      B(j,k) /= m_matrix(k,k);
    }

    //  Solve transpose(L) * X = Y
    for (int k = n-1; k >= j; k--) {
      for (int i = k+1; i < n; i++) 
        B(j,k) -= B(j,i) * m_matrix(i,k);

      B(j,k) /= m_matrix(k,k);
    }
  }

  B.copyUpperTri();
  return B;
}

template<typename MatrixType>
inline MatrixType Cholesky<MatrixType>::solve(const MatrixType& b) const
{
  MatrixType x(b);
  int n = m_matrix.rows();

  //  Solve L*y = b
  for (int k = 0; k < n; k++) {
    for (int i = 0; i < k; i++) 
      x(k) -= x(i) * m_matrix(k,i);

    x(k) /= m_matrix(k,k);
  }

  //  Solve transpose(L) * X = Y
  for (int k = n-1; k >= 0; k--) {
    for (int i = k+1; i < n; i++) 
      x(k) -= x(i) * m_matrix(i,k);

    x(k) /= m_matrix(k,k);
  }
  return b;
}

#endif
