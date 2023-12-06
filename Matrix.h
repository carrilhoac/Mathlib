
#ifndef _Math_Matrix_h
#define _Math_Matrix_h 1

template<typename MatrixLayout>
class MatrixGeneric {
public:
  typedef typename MatrixLayout::DataType DataType;

protected:
  MatrixLayout m_matrix;

public:
  inline MatrixGeneric& resize(int numRows, int numCols) 
  { m_matrix.resize(numRows, numCols);  return *this; }

  inline int rows() const { return m_matrix.rows(); }
  inline int cols() const { return m_matrix.cols(); }
  inline const DataType& operator () (int i, int j) const { return m_matrix(i,j); }
  inline       DataType& operator () (int i, int j)       { return m_matrix(i,j); }
  inline const DataType& operator () (int i) const        { return m_matrix(i); }
  inline       DataType& operator () (int i)              { return m_matrix(i); }

public:
  // For dynamic sized matrices it defaults to a 1 x 1 (or numRows x 1 {vector} if only the numRows is passed)
  // For fixed sized matrices it defaults to the fixed size (4 x 4, for instance)
  inline MatrixGeneric(
    int numRows = MatrixLayout::rowsDefault, 
    int numCols = MatrixLayout::colsDefault) 
  { resize(numRows, numCols); }

  template<typename OtherLayout>
  inline MatrixGeneric(const MatrixGeneric<OtherLayout>& B)
  { copy(B); }

  template<typename OtherLayout>
  inline MatrixGeneric& operator = (const MatrixGeneric<OtherLayout>& B) 
  { return copy(B); }

  inline MatrixGeneric(const std::initializer_list<std::initializer_list<DataType>>& B)
  { setFromInitList(B); }

  inline MatrixGeneric& operator = (
    const std::initializer_list<std::initializer_list<DataType>>& B) 
  { return setFromInitList(B); }

  MatrixGeneric operator - () const;

  MatrixGeneric t() const 
  { return transpose(*this); }

  MatrixGeneric& setZero();
  MatrixGeneric& setIdentity();

  MatrixGeneric getRow(int indexRow) const;
  MatrixGeneric getCol(int indexCol) const;

  // Warning: these methods will fail for sparse matrices
  MatrixGeneric& swapRows(int row1, int row2);
  MatrixGeneric& swapCols(int col1, int col2);

  MatrixGeneric& clearLowerTri();
  MatrixGeneric& clearUpperTri();

  MatrixGeneric& copyUpperTri();
  MatrixGeneric& copyLowerTri();

  MatrixGeneric transposeMultiplySelf() const;

  bool isSquare   () const;
  bool isZero     (DataType tol = eps<DataType>()) const;
  bool isIdentity (DataType tol = eps<DataType>()) const;
  bool isSymmetric(DataType tol = eps<DataType>()) const;
  bool isDiagonal (DataType tol = eps<DataType>()) const;
  bool isTriUpper (DataType tol = eps<DataType>()) const;
  bool isTriLower (DataType tol = eps<DataType>()) const;
  bool isOrthogonal (DataType tol = eps<DataType>()) const;
  //bool isOrthonormal(DataType tol = eps<DataType>()) const;

private:
  template<typename OtherLayout>
  inline MatrixGeneric& copy(const MatrixGeneric<OtherLayout>& B) 
  {
    int m = B.rows();
    int n = B.cols();

    resize(m, n);

    for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      (*this)(i,j) = DataType(B(i,j));

    return *this;
  }

  inline MatrixGeneric& setFromInitList(
    const std::initializer_list<std::initializer_list<DataType>>& B) 
  {
    int nRows = int(B.size());
    int nCols = int(B.begin()->size());
    resize( nRows, nCols );

    int i(0), j(0);

    // this makes me sad too :(
    // TODO: change this to use iterators
    for (const auto& matRow : B) {
    for (const auto& matEle : matRow) 
      { 
        (*this)(i,j) = matEle;  
        j++; 
      }
      j=0; 
      i++;
    }
    return *this;
  }
};

template<typename MatrixType>
inline MatrixGeneric<MatrixType> 
MatrixGeneric<MatrixType>::operator - () const 
{
  const int m = rows();
  const int n = cols();

  MatrixGeneric B(m, n);

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    B(i,j) = -(*this)(i,j);

  return B;
}

template<typename MatrixType>
inline MatrixGeneric<MatrixType>&
MatrixGeneric<MatrixType>::setZero()
{
  const int m = rows();
  const int n = cols();

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    (*this)(i,j) = DataType(0);

  return *this;
}

template<typename MatrixType>
inline MatrixGeneric<MatrixType>&
MatrixGeneric<MatrixType>::setIdentity()
{
  const int m = rows();
  const int n = cols();

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    (*this)(i,j) = i == j ? DataType(1) : DataType(0);

  return *this;
}

#if 0
template<typename MatrixType>
inline MatrixGeneric<MatrixType> 
MatrixGeneric<MatrixType>::transpose() const 
{
  const int m = rows();
  const int n = cols();

  MatrixGeneric<MatrixType> T(n, m);

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    T(j,i) = (*this)(i,j);

  return T;
}
#endif

template<typename MatrixType>
inline MatrixGeneric<MatrixType>& 
MatrixGeneric<MatrixType>::swapRows(int row1, int row2) 
{
  const int n = cols();
  for (int j = 0; j < n; j++)
    std::swap( (*this)(row1,j), (*this)(row2,j) );
    
  return *this;
}
template<typename MatrixType>
inline MatrixGeneric<MatrixType>& 
MatrixGeneric<MatrixType>::swapCols(int col1, int col2) 
{
  const int m = rows();
  for (int i = 0; i < m; i++)
    std::swap( (*this)(i,col1), (*this)(i,col2) );

  return *this;
}

template<typename MatrixType>
inline MatrixGeneric<MatrixType> 
MatrixGeneric<MatrixType>::getRow(int indexRow) const 
{
  const int m = rows();
  const int n = cols();

  MatrixGeneric<MatrixType> R(1, n);

  // Python influence 
  if (indexRow < 0) 
    indexRow += m;

  for (int j = 0; j < n; j++)
    R(0,j) = (*this)(indexRow, j);

  return R;
}

template<typename MatrixType>
inline MatrixGeneric<MatrixType> 
MatrixGeneric<MatrixType>::getCol(int indexCol) const 
{
  const int m = rows();
  const int n = cols();

  MatrixGeneric<MatrixType> C(m, 1);

  // Python influence 
  if (indexCol < 0) 
    indexCol += n;

  for (int i = 0; i < m; i++)
    C(i,0) = (*this)(i, indexCol);

  return C;
}

template<typename MatrixType>
MatrixGeneric<MatrixType>& 
MatrixGeneric<MatrixType>::clearLowerTri() 
{
  const int n = rows();

  for (int i = 0; i < n; i++)
  for (int j = 0; j < i; j++)
    (*this)(i,j) = DataType(0);

  return *this;
}
template<typename MatrixType>
MatrixGeneric<MatrixType>& 
MatrixGeneric<MatrixType>::clearUpperTri() 
{
  const int n = rows();

  for (int i =   0; i < n; i++)
  for (int j = i+1; j < n; j++)
    (*this)(i,j) = DataType(0);

  return *this;
}
template<typename MatrixType>
MatrixGeneric<MatrixType>& 
MatrixGeneric<MatrixType>::copyUpperTri() 
{
  const int n = rows();

  for (int i = 0; i < n; i++)
  for (int j = 0; j < i; j++)
    (*this)(i,j) = (*this)(j,i);

  return *this;
}
template<typename MatrixType>
MatrixGeneric<MatrixType>& 
MatrixGeneric<MatrixType>::copyLowerTri() 
{
  const int n = rows();

  for (int i =   0; i < n; i++)
  for (int j = i+1; j < n; j++)
    (*this)(i,j) = (*this)(j,i);

  return *this;
}

template<typename MatrixType>
inline MatrixGeneric<MatrixType>
MatrixGeneric<MatrixType>::transposeMultiplySelf() const 
{
  const int m = rows();
  const int n = cols();

  MatrixGeneric<MatrixType> B(n, n);

  for (int k = 0; k <  m; k++)
  for (int i = 0; i <  n; i++)
  for (int j = 0; j <= i; j++)
    B(i,j) += (*this)(k,i) * (*this)(k,j);

  B.copyLowerTri();
  return B;
}

template<typename MatrixType>
inline bool 
MatrixGeneric<MatrixType>::isSquare() const 
{
  return rows() == cols();
}

template<typename MatrixType>
inline bool 
MatrixGeneric<MatrixType>::isZero(typename MatrixType::DataType tol) const 
{
  const int m = rows();
  const int n = cols();

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    if (std::fabs((*this)(i,j)) > tol)
      return false;

  return true;
}

template<typename MatrixType>
inline bool 
MatrixGeneric<MatrixType>::isIdentity(typename MatrixType::DataType tol) const 
{
  const int m = rows();
  const int n = cols();

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
  {
    if (i != j && (std::fabs((*this)(i,j)) > tol))
      return false;
    if (i == j && (std::fabs((*this)(i,j) - DataType(1) ) > tol))
      return false;
  }

  return true;
}

template<typename MatrixType>
inline bool 
MatrixGeneric<MatrixType>::isSymmetric(typename MatrixType::DataType tol) const 
{
  const int m = rows();
// const int n = cols();

  for (int i = 0; i < m; i++)
  for (int j = 0; j < i; j++)
    if (std::fabs((*this)(i,j) - (*this)(j,i)) > tol)
      return false;

  return true;
}

template<typename MatrixType>
inline bool 
MatrixGeneric<MatrixType>::isDiagonal(typename MatrixType::DataType tol) const
{
  if (! isSquare())
    return false;

  const int m = rows();
  const int n = cols();

  for (int i = 0; i < m; i++)
  for (int j = 0; j < n; j++)
    if (i != j && std::fabs((*this)(i,j)) > tol)
      return false;

  return true;
}

template<typename MatrixType>
inline bool 
MatrixGeneric<MatrixType>::isTriUpper(typename MatrixType::DataType tol) const
{
  if (! isSquare())
    return false;

  const int m = rows();
  const int n = cols();

  for (int i = 0  ; i < m; i++)
  for (int j = i+1; j < n; j++)
    if (std::fabs((*this)(i,j)) > tol)
      return false;

  return true;
}

template<typename MatrixType>
inline bool 
MatrixGeneric<MatrixType>::isTriLower(typename MatrixType::DataType tol) const
{
  if (! isSquare())
    return false;

  const int m = rows();
  const int n = cols();

  for (int i = 0; i < m; i++)
  for (int j = 0; j < i; j++)
    if (std::fabs((*this)(i,j)) > tol)
      return false;

  return true;
}

template<typename MatrixType>
inline bool 
MatrixGeneric<MatrixType>::isOrthogonal(typename MatrixType::DataType tol) const
{
  if (! isSquare())
    return false;

  return transposeMultiplySelf().isIdentity(tol);
}

#endif


