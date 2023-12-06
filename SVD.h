
#ifndef _Math_SVD_h
#define _Math_SVD_h 1

// Singular Value Decomposition
//
//   For an m-by-n matrix A with m >= n, the singular value decomposition is
//   an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
//   an n-by-n orthogonal matrix V so that A = U*S*V'.
//
// The singular values, sigma[k], are ordered so that
//   sigma[0] >= sigma[1] >= ... >= sigma[n-1].
//
// The singular value decompostion always exists, so the constructor will
//   never fail.  The matrix condition number and the effective numerical
//   rank can be computed from this decomposition.

template<typename MatrixType>
class SVD {
protected:
  typedef typename MatrixType::DataType DataType;
  typedef typename std::vector<DataType> VectorType;

protected:
  MatrixType m_matU;
  MatrixType m_matV;
  VectorType m_sigma;

public:
  SVD(const MatrixType& A)
  { factorize(A); }

  int rank(DataType tol = eps<DataType>()) const;

  MatrixType getU() const { return m_matU; }
  MatrixType getV() const { return m_matV; }
  MatrixType getS() const;

  DataType cond() const 
  { return m_sigma.front() / m_sigma.back(); }

  DataType norm2() const 
  { return m_sigma.front(); }
  
  MatrixType nullspace() const  // todo probably check if last singular value is null as well  ?
  { return m_matV.getCol(-1); } // last column of matrix V

private:
  void factorize(const MatrixType& A);
};

template<typename MatrixType>
int SVD<MatrixType>::rank(DataType tol) const 
{ 
  int rank_ = 0;
  int n = int(m_sigma.size());
  for (int i = 0; i < n; i++) {
    if (m_sigma[i] > tol) {
      rank_++;
    }
  }
  return rank_; 
}

template<typename MatrixType>
MatrixType SVD<MatrixType>::getS() const 
{ 
  int n = int(m_sigma.size());

  MatrixType S(n,n);
  S.setZero();

  for (int i = 0; i < n; i++)
    S(i,i) = m_sigma[i];

  return S;
}

template<typename MatrixType>
void SVD<MatrixType>::factorize(const MatrixType& A)
{
  const int m = A.rows();
  const int n = A.cols();
  const int nu = std::min(m,n);

  const bool wantu = true;
  const bool wantv = true;

  m_matU.resize(m,nu);
  m_matU.setZero();
  m_sigma.resize(nu);
  m_matV.resize(n,n);
  
  VectorType e(n);
  VectorType work(m);
  MatrixType matA(A);

  int i=0, j=0, k=0;

  // Reduce A to bidiagonal form, storing the diagonal elements
  // in s and the super-diagonal elements in e.
  const int nct = std::min(m-1, n);
  const int nrt = std::max(0, std::min(n-2,m));

  for (k = 0; k < std::max(nct, nrt); k++) {
    if (k < nct) {

      // Compute the transformation for the k-th column and
      // place the k-th diagonal in m_sigma[k].
      // Compute 2-norm of k-th column without under/overflow.
      m_sigma[k] = DataType(0);
      for (i = k; i < m; i++) {
        m_sigma[k] = std::hypot(m_sigma[k], matA(i,k));
      }
      if (m_sigma[k] != DataType(0)) { 
        if (matA(k,k) < DataType(0)) {
          m_sigma[k] = -m_sigma[k];
        }
        for (i = k; i < m; i++) {
          matA(i,k) /= m_sigma[k];
        }
        matA(k,k) += DataType(1);
      }
      m_sigma[k] = -m_sigma[k];
    }

    for (j = k+1; j < n; j++) {
      if ((k < nct) && (m_sigma[k] != DataType(0))) {

        // Apply the transformation.
        DataType t(0);
        for (i = k; i < m; i++) {
          t += matA(i,k) * matA(i,j);
        }
        t = -t / matA(k,k);
        for (i = k; i < m; i++) {
          matA(i,j) += t * matA(i,k);        
        }
      }

      // Place the k-th row of A into e for the
      // subsequent calculation of the row transformation.
      e[j] = matA(k,j);
    }

    // Place the transformation in U for subsequent back multiplication.
    if (wantu & (k < nct)) {
      for (i = k; i < m; i++) {
        m_matU(i,k) = matA(i,k);
      }
    }

    if (k < nrt) {
      // Compute the k-th row transformation and place the
      // k-th super-diagonal in e[k].
      // Compute 2-norm without under/overflow.
      e[k] = DataType(0);
      for (i = k+1; i < n; i++) {
        e[k] = std::hypot(e[k], e[i]);
      }
      if (e[k] != DataType(0)) {
        if (e[k+1] < DataType(0)) {
          e[k] = -e[k];
        }        
        for (i = k+1; i < n; i++) {
          e[i] /= e[k];
        }
        e[k+1] += DataType(1);
      }

      e[k] = -e[k];
      if ((k+1 < m) & (e[k] != DataType(0))) {
        // Apply the transformation.
        for (i = k+1; i < m; i++) {
          work[i] = DataType(0);
        }
        
        for (j = k+1; j < n; j++) {
          for (i = k+1; i < m; i++) {
            work[i] += e[j] * matA(i,j);
          }
        }

        for (j = k+1; j < n; j++) {
          DataType t = -e[j] / e[k+1];
          
          for (i = k+1; i < m; i++)  {
            matA(i,j) += t * work[i];
          }
        }
      }

      // Place the transformation in V for subsequent back multiplication.
      if (wantv) {
        for (i = k+1; i < n; i++) {
          m_matV(i,k) = e[i];
        }
      }
    }
  }

  // Set up the final bidiagonal matrix or order p.
  int p = std::min(n, m+1);
  if (nct < n) {
    m_sigma[nct] = matA(nct,nct);
  }
  if (m < p)  {
    m_sigma[p-1] = DataType(0);
  }
  if ((nrt+1) < p) {
    e[nrt] = matA(nrt,p-1);
  }
  e[p-1] = DataType(0);

  // If required, generate U.
  if (wantu) {
    for (j = nct; j < nu; j++) {
      for (i = 0; i < m; i++) {
        m_matU(i,j) = DataType(0);
      }
      m_matU(j,j) = DataType(1);
    }
    for (k = nct-1; k >= 0; k--) {
      if (m_sigma[k] != DataType(0)) {
        for (j = k+1; j < nu; j++) {
          DataType t(0);
          for (i = k; i < m; i++) {
            t += m_matU(i,k) * m_matU(i,j); 
          }
          t = -t / m_matU(k,k);

          for (i = k; i < m; i++) {
            m_matU(i,j) += t * m_matU(i,k);
          }
        }
        for (i = k; i < m; i++ ) {
          m_matU(i,k) = -m_matU(i,k);
        }
        
        m_matU(k,k) += DataType(1);

        for (i = 0; i < k-1; i++) {
          m_matU(i,k) = DataType(0);
        }
      } 
      else {
        for (i = 0; i < m; i++) {
          m_matU(i,k) = DataType(0);
        }
        
        m_matU(k,k) = DataType(1);
      }
    }
  }

  // If required, generate V.
  if (wantv) {
    for (k = n-1; k >= 0; k--) {
      if ((k < nrt) & (e[k] != DataType(0))) {
        for (j = k+1; j < nu; j++) {
          DataType t(0);

          for (i = k+1; i < n; i++) {
            t += m_matV(i,k) * m_matV(i,j); 
          }
          
          t = -t / m_matV(k+1,k);

          for (i = k+1; i < n; i++) {
            m_matV(i,j) += t * m_matV(i,k);
          }
        }
      }
      for (i = 0; i < n; i++) {
        m_matV(i,k) = DataType(0);
      }
      
      m_matV(k,k) = DataType(1);
    }
  }

  // Main iteration loop for the singular values.
  int pp = p-1;
  int iter = 0;
  
  //DataType tol = std::pow(2.0, -52.0);
  DataType tol = eps<DataType>();

  while (p > 0) {
    //int k = 0;
    int k = 0;
    int kase = 0;

    // Here is where a test for too many iterations would go.

    // This section of the program inspects for
    // negligible elements in the s and e arrays.  On
    // completion the variables kase and k are set as follows.

    // kase = 1     if s(p) and e[k-1] are negligible and k<p
    // kase = 2     if s(k) is negligible and k<p
    // kase = 3     if e[k-1] is negligible, k<p, and
    //              s(k), ..., s(p) are not negligible (QR step).
    // kase = 4     if e(p-1) is negligible (convergence).

    for (k = p-2; k >= -1; k--) {
      if (k == -1) {
        break;
      }
      
      if (std::fabs(e[k]) <= tol*(std::fabs(m_sigma[k]) + std::fabs(m_sigma[k+1]))) {
        e[k] = DataType(0);
        break;
      }
    }

    if (k == (p-2)) {
      kase = 4;
    }
    else {
      int ks;
      for (ks = p-1; ks >= k; ks--) {
        if (ks == k) {
          break;
        }
        
        DataType t = 
          (ks != p   ? std::fabs(e[ks  ]) : DataType(0)) + 
          (ks != k+1 ? std::fabs(e[ks-1]) : DataType(0));

        if (std::fabs(m_sigma[ks]) <= tol*t) {
          m_sigma[ks] = DataType(0);
          break;
        }
      }

      if (ks == k) {
        kase = 3;
      }
      else if (ks == (p-1)) {
        kase = 1;
      }
      else {
        kase = 2;
        k = ks;
      }
    }
    k++;

    // Perform the task indicated by kase.
    switch (kase) {

      // Deflate negligible s(p).
      case 1: {
        DataType f = e[p-2];
        e[p-2] = DataType(0);
        for (j = p-2; j >= k; j--) {
          DataType t(std::hypot(m_sigma[j], f));
          DataType cs(m_sigma[j] / t);
          DataType sn(f/t);
          m_sigma[j] = t;
          if (j != k) {
            f      = -sn * e[j-1];
            e[j-1] =  cs * e[j-1];
          }
          if (wantv) {
            for (i = 0; i < n; i++) {
              t = cs*m_matV(i,j) + sn*m_matV(i,p-1);
              m_matV(i,p-1) = -sn*m_matV(i,j) + cs*m_matV(i,p-1);
              m_matV(i,j) = t;
            }
          }
        }
      }
      break;

      // Split at negligible s(k).
      case 2: {
        DataType f = e[k-1];
        e[k-1] = DataType(0);
        for (j = k; j < p; j++) {
          DataType t(std::hypot(m_sigma[j],f));
          DataType cs(m_sigma[j]/t);
          DataType sn(f/t);
          m_sigma[j] = t;
          f    = -sn*e[j];
          e[j] =  cs*e[j];
          if (wantu) {
            for (i = 0; i < m; i++) {
              t = cs*m_matU(i,j) + sn*m_matU(i,k-1);
              m_matU(i,k-1) = -sn*m_matU(i,j) + cs*m_matU(i,k-1);
              m_matU(i,j) = t;
            }
          }
        }
      }
      break;

      // Perform one QR step.
      case 3: {

        // Calculate the shift.
        DataType scale = 
          std::max(   // just amazing, don't you think?
            std::max(
              std::max(
                std::max(
                  std::fabs(m_sigma[p-1]),
                  std::fabs(m_sigma[p-2])),
                std::fabs(e[p-2])),
              std::fabs(m_sigma[k])),
            std::fabs(e[k]));

        DataType sp   = m_sigma[p-1]/scale;
        DataType spm1 = m_sigma[p-2]/scale;
        DataType epm1 = e[p-2]/scale;
        DataType sk   = m_sigma[k]/scale;
        DataType ek   = e[k]/scale;
        DataType b    = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/DataType(2);
        DataType c    = (sp*epm1)*(sp*epm1);
        DataType shift(0);
        if ((b != DataType(0)) || (c != DataType(0))) {
          shift = std::sqrt(b*b + c);
          if (b < DataType(0)) {
            shift = -shift;
          }
          shift = c/(b + shift);
        }
        DataType f = (sk + sp)*(sk - sp) + shift;
        DataType g = sk*ek;

        // Chase zeros.
        for (j = k; j < p-1; ++j) {
          DataType t = std::hypot(f,g);
          DataType cs = f/t;
          DataType sn = g/t;
          if (j != k) {
            e[j-1] = t;
          }
          f = cs*m_sigma[j] + sn*e[j];
          e[j] = cs*e[j] - sn*m_sigma[j];
          g = sn*m_sigma[j+1];
          m_sigma[j+1] = cs*m_sigma[j+1];
          if (wantv) {
            for (i = 0; i < n; ++i) {
              t = cs*m_matV(i,j) + sn*m_matV(i,j+1);
              m_matV(i,j+1) = -sn*m_matV(i,j) + cs*m_matV(i,j+1);
              m_matV(i,j) = t;
            }
          }
          t = std::hypot(f,g);
          cs = f/t;
          sn = g/t;
          m_sigma[j] = t;
          f = cs*e[j] + sn*m_sigma[j+1];
          m_sigma[j+1] = -sn*e[j] + cs*m_sigma[j+1];
          g = sn*e[j+1];
          e[j+1] = cs*e[j+1];
          if (wantu && (j < m-1)) {
            for (i = 0; i < m; ++i) {
              t = cs*m_matU(i,j) + sn*m_matU(i,j+1);
              m_matU(i,j+1) = -sn*m_matU(i,j) + cs*m_matU(i,j+1);
              m_matU(i,j) = t;
            }
          }
        }
        e[p-2] = f;
        iter = iter + 1;
      }
      break;
      
      // Convergence.
      case 4: {
        
        // Make the singular values positive.
        if (m_sigma[k] <= DataType(0)) {
          m_sigma[k] = (m_sigma[k] < DataType(0) ? -m_sigma[k] : DataType(0));
          if (wantv) {
            for (i = 0; i <= pp; i++) {
              m_matV(i,k) = -m_matV(i,k);
            }
          }
        }
   
        // Order the singular values.
        while (k < pp) {
          if (m_sigma[k] >= m_sigma[k+1]) {
            break;
          }
          
          std::swap( m_sigma[k], m_sigma[k+1] );

          if (wantv && (k < n-1)) {
            m_matV.swapCols(k, k+1);
          }

          if (wantu && (k < m-1)) {
            m_matU.swapCols(k, k+1);
          }
          k++;
        }
        iter = 0;
        p--;
      }
      break;
    } 
  }
}

#endif

