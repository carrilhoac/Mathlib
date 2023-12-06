
#ifndef _Math_QR_h
#define _Math_QR_h 1

// Householder QR decompisition A = Q * R where:
//   A  m-by-n  (m >= n) matrix
//   Q  m-by-n  orthogonal matrix
//   R  n-by-n  upper triangular matrix
// The decompostion always exists, even if the matrix does not have
// full rank, so the constructor will never fail. 
template<typename MatrixType>
class QR {
protected:
  typedef typename MatrixType::DataType DataType;
  typedef typename std::vector<DataType> VectorType;

protected:
  MatrixType m_matrix;
  VectorType m_Rdiag;
  bool m_isFullRank;

public:
  QR(const MatrixType& A)
    : m_matrix(A)
    , m_Rdiag(A.cols())
    , m_isFullRank(false)
  { factorize(); }

  MatrixType getQ() const;
  MatrixType getR() const;

  MatrixType solve(const MatrixType& b) const;

private:
  void factorize();
};

template<typename MatrixType>
void QR<MatrixType>::factorize() 
{
  int m = m_matrix.rows();
  int n = m_matrix.cols();

  for (int k = 0; k < n; k++) {
    // compute 2-norm of k-th column 
    DataType nrm(0);

    for (int i = k; i < m; i++)
      nrm = std::hypot(nrm, m_matrix(i,k));

    if (nrm != DataType(0)) {
      // k-th Householder vector
      if (m_matrix(k,k) < DataType(0))
        nrm = -nrm;

      for (int i = k; i < m; i++)
        m_matrix(i,k) /= nrm;

      m_matrix(k,k) += DataType(1);

      // apply transformation to remaining columns
      for (int j = k+1; j < n; j++) {
        DataType s(0);

        for (int i = k; i < m; i++)
          s += m_matrix(i,k) * m_matrix(i,j);
        
        s = -s / m_matrix(k,k);

        for (int i = k; i < m; i++)
          m_matrix(i,j) += s * m_matrix(i,k);
      }
    }
    m_Rdiag[k] = -nrm;
  }

  // TODO: investigate this... maybe should use a epsilon(?)
  for (int j = 0; j < n; j++)
    if (m_Rdiag[j] == DataType(0))
      m_isFullRank = false;
}

template<typename MatrixType>
MatrixType QR<MatrixType>::getQ() const 
{
  int m = m_matrix.rows();
  int n = m_matrix.cols();
  MatrixType Q(m,n);

  for (int k = n-1; k >= 0; k--) {
    Q(k,k) = DataType(1);

    for (int j = k; j < n; j++) {
      if (m_matrix(k,k) != DataType(0)) {
        DataType s(0);

        for (int i = k; i < m; i++) 
          s += m_matrix(i,k) * Q(i,j);

        s = -s / m_matrix(k,k);

        for (int i = k; i < m; i++) 
          Q(i,j) += s * m_matrix(i,k);
      }
    }
  }
  return Q;
}

template<typename MatrixType>
MatrixType QR<MatrixType>::getR() const 
{
  int n = m_matrix.cols();
  MatrixType R(n,n);

  for (int i = 0; i < n; i++)
  for (int j = i; j < n; j++)
    R(i,j) = i == j ? m_Rdiag[i] : m_matrix(i,j);

  return R;
}

// Least squares solution of A*x = b
//  b  is m-by-1  vector right hand side.
//  x  is n-by-1  vector that minimizes the 2-norm of Q*R*X-B.
template<typename MatrixType>
MatrixType QR<MatrixType>::solve(const MatrixType& b) const 
{
  int m = m_matrix.rows();
  int n = m_matrix.cols();

  //if (! m_isFullRank)
  //  return MatrixType(n);

  MatrixType x_ = b;
  MatrixType x(n);

  // Compute Y = transpose(Q)*b
  for (int k = 0; k < n; k++) {
    DataType s(0); 
    for (int i = k; i < m; i++) {
      s += m_matrix(i,k) * x_(i);
    }
    s = -s / m_matrix(k,k);
    for (int i = k; i < m; i++) {
      x_(i) += s*m_matrix(i,k);
    }
  }
  // Solve R*X = Y;
  for (int k = n-1; k >= 0; k--) {
    x_(k) /= m_Rdiag[k];
    for (int i = 0; i < k; i++) {
      x_(i) -= x_(k) * m_matrix(i,k);
    }
  }

  // return n x 1 portion of X
  for (int i = 0; i < n; i++)
    x(i) = x_(i);

  return x_;
}

#endif

