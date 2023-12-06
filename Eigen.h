
#ifndef _Math_Eigen_h
#define _Math_Eigen_h 1

//  Computes eigenvalues and eigenvectors of a real (non-complex)
//  matrix. Matrix must be square !
//   
//  If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
//  diagonal and the eigenvector matrix V is orthogonal. That is,
//	the diagonal values of D are the eigenvalues, and
//  V*V' = I, where I is the identity matrix. The columns of V 
//  represent the eigenvectors in the sense that A*V = V*D.
//
//  Note: If A is symmetric and positive definite, this is equal to 
//  the SVD decomposition (Spectral theorem)
//
//  If A is non-symmetric, then A = V*D*inv(V) 

template<typename MatrixType>
class Eigen {
protected:
  typedef typename MatrixType::DataType DataType;
  typedef typename std::vector<DataType> VectorType;

protected:
  int n;
  VectorType d;   // eigenvalues  <real> part
  VectorType e;   // eigenvalues  <img>  part
  MatrixType V;   // eigenvectors (columns of V)

public:
  MatrixType getD() const;
  MatrixType getV() const 
  { return V; }

  Eigen(const MatrixType& A) 
    : n (A.rows())
    , d (n)
    , e (n)
    , V (n, n)
  {
    if (A.isSymmetric())
    {
      V = A;
      tred2(); // Symmetric Householder reduction to tridiagonal form
      tql2();  // Diagonalize (Symmetric tridiagonal QL algorithm)
    } 
    else 
    {
      // Matrix must be square
      MatrixType H = A;    // non-symmetric Hessenberg form
      VectorType ort(n);   // non-symmetric help storage
      
      orthes(H, ort); // Reduce to Hessenberg form.
      hqr2(H, ort);   // Reduce Hessenberg to real Schur form.
    }

    eigsort();
  } 

private:
  void eigsort(); 

  void tred2();
  void tql2();
  
  void orthes(MatrixType& H, VectorType& ort);
  void hqr2(MatrixType& H, VectorType& ort);
};

template<typename MatrixType>
MatrixType Eigen<MatrixType>::getD() const 
{
  MatrixType D(n,n);
  D.setZero();

  for (int i = 0; i < n; i++) {
    D(i,i) = d[i];

    if (e[i] > DataType(0)) {
      D(i,i+1) = e[i];
    }
    if (e[i] < DataType(0)) {
      D(i,i-1) = e[i];
    }
  }

  return D;
}

template<typename MatrixType>
void Eigen<MatrixType>::eigsort() 
{
  // Init permutation vector
  std::vector<int> perm(n);
  for (int i = 0; i < n; i++)
    perm[i] = i;

  // Sort permutation vector based on eigen value's <real> part
  // from larger to smaller  
  std::sort( perm.begin(), perm.end(),
    [&](const int& a, const int& b) 
    { return d[a] > d[b]; } );
  
  // Make permutations
  VectorType d_(n);
  VectorType e_(n);
  MatrixType V_(n,n);

  // Doing it like this duplicates the memory footprint for a little while
  // however, it minimizes the number of 'movs', which I think is better
  for (int j = 0; j < n; j++) {
    const int k = perm[j];
    d_[j] = d[k];
    e_[j] = e[k];

    for (int i = 0; i < n; i++) 
      V_(i,j) = V(i,k);
  }
  d = d_;
  e = e_;
  V = V_;
}

#endif 

