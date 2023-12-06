
#ifndef _Math_LU_h
#define _Math_LU_h 1

// Partial pivoting LU decomposition (only rows are pivoted, not columns)
//    P*A = L*U
//
//  P*A = L*U   =>  (P^-1)*P*A = (P^-1)*P*L*U   =>   A = (P^-1)*P*L*U
//  Since P is orthonormal:   P^-1 = P^T    and    A = (P^T)*L*U
//
// For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
//  unit lower triangular matrix L, an n-by-n upper triangular matrix U,
//  and a permitation matrix P (here stored as a vector), but it can 
//  be retrieved by getP()
//
template<typename MatrixType>
class LU {
protected:
  typedef typename MatrixType::DataType DataType;
  typedef typename std::vector<DataType> VectorType;
  typedef typename std::vector<int> PivotType;

protected:
public:
  MatrixType m_matrix;
  PivotType m_pivot;
  DataType m_det;
  
public:
  LU(const MatrixType& A) 
    : m_matrix(A)
    , m_pivot(A.rows())
    , m_det(0)
  {
    factorize();
  }
 
  MatrixType getP() const; // permutation matrix
  MatrixType getL() const;
  MatrixType getU() const;

  DataType det() const { return m_det; }
  MatrixType solve(const MatrixType& b) const;
  MatrixType inverse() const;
  
private:
  void factorize();
};

template<typename MatrixType>
MatrixType LU<MatrixType>::getP() const 
{
  int n = int(m_pivot.size());

  MatrixType P(n,n);
  P.setZero();

  for (int i = 0; i < n; i++)
    P(i, m_pivot[i]) = DataType(1);

  return P;
}

template<typename MatrixType>
MatrixType LU<MatrixType>::getL() const 
{
  int m = m_matrix.rows();
  int n = m_matrix.cols();

  MatrixType L(m,n);
  L.setIdentity();

  for (int i = 0; i < m; i++) { 
  for (int j = 0; j < i; j++) { 
    L(i,j) = m_matrix(i,j);
  } }

  return L;
}

template<typename MatrixType>
MatrixType LU<MatrixType>::getU() const 
{
  int n = m_matrix.cols();

  MatrixType U(n,n);
  U.setZero();

  for (int i = 0; i < n; i++) { 
  for (int j = i; j < n; j++) { 
    U(i,j) = m_matrix(i,j);
  } }

  return U;
}

template<typename MatrixType>
void LU<MatrixType>::factorize()
{ 
  int m = m_matrix.rows();
  int n = m_matrix.cols();
  int q = std::min(m,n);

  int i, j, k, p;
  DataType s(0), t(0), d(0), w(1);

  // initialize array for partial pivoting
  for (i = 0; i < m; i++)
    m_pivot[i] = i;
  
  for (i = 0; i < q; i++) {

    // find pivot
    p = i;
    s = std::fabs( m_matrix(i,i) );  
    for (j = i+1; j < m; j++) { 
      t = std::fabs( m_matrix(j,i) );
      if (t > s) 
        { p = j; s = t; }
    }

    // swap if necessary
    if (p != i) {
      w = -w; // swap det sign
      std::swap(m_pivot[i], m_pivot[p]);
      m_matrix.swapRows(i, p);
    }

    // Gaussian elimination
    if (i < m) {
      d = DataType(1) / m_matrix(i,i);
      for (j = i+1; j < m; j++) 
        m_matrix(j,i) *= d;        
    }
    if (i < (q-1)) {
      for (j = i+1; j < m; j++){
        d = m_matrix(j,i);

        for (k = i+1; k < n; k++) 
          m_matrix(j,k) -= d * m_matrix(i,k);
      }
    }
  }

  // computing the determinant
  for (i = 0; i < n; i++) {
    w *= m_matrix(i,i);
  }
  m_det = w;
}

template<typename MatrixType>
MatrixType LU<MatrixType>::solve(const MatrixType& b) const 
{
  int n = m_matrix.cols();
  MatrixType x(n);

  for (int i = 0; i < n; i++) {
    DataType curB = b[ m_pivot[i] ];

    for (int j = 0; j < i; j++) 
      curB -= m_matrix(i,j) * x[j];
    
    x[i] = curB;
  }

  for (int i = n-1; i >= 0; i--) {
    DataType curX = x[i];

    for (int j = i+1; j < n; j++)
      curX -= m_matrix(i,j) * x[j];
    
    x[i] = curX / m_matrix(i,i);
  }
}

template<typename MatrixType>
MatrixType LU<MatrixType>::inverse() const 
{
//  Similar strategy to Cholesky: Here we solve the inverse of A, 
//  but stored as inv(A)^T However, as the input matrix is not 
//  ensured to be symmetric, we need to return the transpose of the result.

  int n = m_matrix.cols();

  std::vector<DataType> b(n);
  MatrixType X(n,n);

  // not needed 
  //for (int k = 0; k < n; k++)
  //  b[k] = DataType(0);

  for (int k = 0; k < n; k++) {
    b[k] = DataType(1);

    for (int i = 0; i < n; i++) {
      DataType curB = b[ m_pivot[i] ];

      for (int j = 0; j < i; j++) 
        curB -= m_matrix(i,j) * X(k,j);
      
      X(k,i) = curB;
    }

    for (int i = n-1; i >= 0; i--) {
      DataType curX = X(k,i);

      for (int j = i+1; j < n; j++)
        curX -= m_matrix(i,j) * X(k,j);
      
      X(k,i) = curX / m_matrix(i,i);
    }

    b[k] = DataType(0);
  }

  return X.t();
}

#endif

