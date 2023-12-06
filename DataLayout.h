
#ifndef _Math_DataLayout_h
#define _Math_DataLayout_h 1

// Specific to internal storage, making later stuff agnostic to it.
// This is needed if we want to change the underlying memory layout
// without having impact in the code elsewhere.
// The disadvantage is that the algorithms won't know the best 
// accessing pattern unless we add some kind of type-trait to hint it.
template<typename Real, int rowsDefault_ = 1, int colsDefault_ = 1>
class DataLayout {
public:
  typedef Real DataType;
  typedef Real ScalarType;

  static const int rowsDefault = rowsDefault_; // constexpr?
  static const int colsDefault = colsDefault_;

public:
  virtual void resize(int numRows_, int numCols_) = 0;

  virtual int rows() const = 0;
  virtual int cols() const = 0;
  virtual const Real& operator () (int i, int j) const = 0;
  virtual       Real& operator () (int i, int j) = 0;
  virtual const Real& operator () (int i) const = 0;
  virtual       Real& operator () (int i) = 0;
};

template<typename Real>
class DenseRowMajor : public DataLayout<Real> {
private:
  int numRows;
  int numCols;

  std::vector<Real> rowMajor;

public:
  inline void resize(int numRows_, int numCols_) {
    numRows = numRows_;
    numCols = numCols_;
    rowMajor.resize( numRows * numCols );
  }
  
  inline int rows() const { return numRows; }
  inline int cols() const { return numCols; }
  inline const Real& operator () (int i, int j) const {
    return rowMajor[ i * numCols + j ];
  }
  inline       Real& operator () (int i, int j){
    return rowMajor[ i * numCols + j ];
  }
  inline const Real& operator () (int i) const {
    return rowMajor[ i ];
  }
  inline       Real& operator () (int i) {
    return rowMajor[ i ];
  }
};

template<typename Real>
class DenseColMajor : public DataLayout<Real> {  
private:
  int numRows;
  int numCols;

  std::vector<Real> colMajor;

public:
  inline void resize(int numRows_, int numCols_) {
    numRows = numRows_;
    numCols = numCols_;
    colMajor.resize( numRows * numCols );
  }
  
  inline int rows() const { return numRows; }
  inline int cols() const { return numCols; }
  inline const Real& operator () (int i, int j) const {
    return colMajor[ j * numRows + i ];
  }
  inline       Real& operator () (int i, int j){
    return colMajor[ j * numRows + i ];
  }
  inline const Real& operator () (int i) const {
    return colMajor[ i ];
  }
  inline       Real& operator () (int i) {
    return colMajor[ i ];
  }
};

template<typename Real, int rowsFixed = 4, int colsFixed = 4>
class DenseRowMajorFixed : public DataLayout<Real, rowsFixed, colsFixed> {  
private:
  int numRows;
  int numCols;

  Real colMajor[ rowsFixed * colsFixed ]; // constexpr

public:
  DenseRowMajorFixed()
    : numRows( this->rowsDefault )
    , numCols( this->colsDefault )
  { }
  // we can resize it to smaller than 4x4, 3x3, or whatever
  // memory will stay the same, however
  inline void resize(int numRows_, int numCols_) {
    numRows = numRows_;
    numCols = numCols_;
  }
  
  inline int rows() const { return numRows; }
  inline int cols() const { return numCols; }
  inline const Real& operator () (int i, int j) const {
    return colMajor[ i * numCols + j ];
  }
  inline       Real& operator () (int i, int j){
    return colMajor[ i * numCols + j ];
  }
  inline const Real& operator () (int i) const {
    return colMajor[ i ];
  }
  inline       Real& operator () (int i) {
    return colMajor[ i ];
  }
};
template<typename Real, int rowsFixed = 4, int colsFixed = 4>
class DenseColMajorFixed : public DataLayout<Real, rowsFixed, colsFixed> {  
private:
  int numRows;
  int numCols;

  Real colMajor[ rowsFixed * colsFixed ]; // constexpr

public:
  DenseColMajorFixed()
    : numRows( this->rowsDefault )
    , numCols( this->colsDefault )
  { }
  // we can resize it to smaller than 4x4, 3x3, or whatever
  // memory will stay the same, however
  inline void resize(int numRows_, int numCols_) {
    numRows = numRows_;
    numCols = numCols_;
  }
  
  inline int rows() const { return numRows; }
  inline int cols() const { return numCols; }
  inline const Real& operator () (int i, int j) const {
    return colMajor[ j * numRows + i ];
  }
  inline       Real& operator () (int i, int j){
    return colMajor[ j * numRows + i ];
  }
  inline const Real& operator () (int i) const {
    return colMajor[ i ];
  }
  inline       Real& operator () (int i) {
    return colMajor[ i ];
  }
};

//template<typename Real, int tileSize = 4> // 4x4 tiles
//class DenseTiled : public DataLayout<Real> {  
//};


namespace Sparse {
//  Sparse matrix with constant layout (most zeros won't change places)
//
//  Warning: These classes should not be used as MatrixGeneric template parameter
//  since they violate some of the operations (swapRows, swapCols, etc).
//  They are barebone indexers, and are only here to be used on very rare 
//  ocasions, for very large matrices. Also, you will have to implement 
//  everything from them, even the basic operations.
template<typename Real>
class _TriangularLower : public DataLayout<Real> {
private:
  int n;
  Real x; // dummy for when accessing elements outside triangular
  Real y;
  std::vector<Real> triLower;

public:
  inline void resize(int numRows_, int numCols_) {
    n = numRows_;
    x = Real(0);
    y = Real(0);
    triLower.resize((n * (n + 1)) >> 1);
  }
  inline int rows() const { return n; }
  inline int cols() const { return n; }

  inline const Real& operator () (int i, int j) const {
    if (j > i) {
      return x;
    }

    return triLower[((i * (i + 1)) >> 1) + j];
  }
  inline       Real& operator () (int i, int j){
    if (j > i) {
      y = Real(0);
      return y;
    }
        
    return triLower[((i * (i + 1)) >> 1) + j];
  }

  // hmm... 
  inline const Real& operator () (int i) const {
    return triLower[ i ];
  }
  inline       Real& operator () (int i) {
    return triLower[ i ];
  }
};

template<typename Real>
class _TriangularUpper : public DataLayout<Real> {
private:
  int n;
  Real x; // dummy for when accessing elements outside triangular
  Real y;
  std::vector<Real> triUpper;

public:
  inline void resize(int numRows_, int numCols_) {
    n = numRows_;
    x = Real(0);
    y = Real(0);
    triUpper.resize((n * (n + 1)) >> 1);
  }
  inline int rows() const { return n; }
  inline int cols() const { return n; }

  inline const Real& operator () (int i, int j) const {
    if (i > j) {
      return x;
    }

    const int index = (((n << 1) * i) - (i * i) -i + (j << 1)) >> 1;
    return triUpper[index];
  }
  inline       Real& operator () (int i, int j){
    if (i > j) {
      y = Real(0);
      return y;
    }
        
    const int index = (((n << 1) * i) - (i * i) -i + (j << 1)) >> 1;
    return triUpper[index];
  }

  // hmm... 
  inline const Real& operator () (int i) const {
    return triUpper[ i ];
  }
  inline       Real& operator () (int i) {
    return triUpper[ i ];
  }
};

template<typename Real>
class _SymmetricRowMajor : public DataLayout<Real> {
private:
  // Square n-by-n symmetric matrix with (n/2)(n+1) elements 
  // only the lower triangular plus the diagonal are stored
  int n;
  std::vector<Real> symRowMajor; 
  
  // Note: Dont change the order on the indexing integer division 
  // operations (bit shift right) or it might cause underflow.
  // The bit shift right must be the last operation
public:
  inline void resize(int numRows_, int numCols_) {
    n = numRows_;
    symRowMajor.resize((n * (n + 1)) >> 1);
  }
  inline int rows() const { return n; }
  inline int cols() const { return n; }
  
  inline const Real& operator () (int i, int j) const {
    const int u = i > j ? i : j; // make sure its on lower tri
    const int v = i > j ? j : i;     
    return symRowMajor[((u * (u + 1)) >> 1) + v];
  }
  inline       Real& operator () (int i, int j){
    const int u = i > j ? i : j; // make sure its on lower tri
    const int v = i > j ? j : i;     
    return symRowMajor[((u * (u + 1)) >> 1) + v];
  }

  // hmm... 
  inline const Real& operator () (int i) const {
    return symRowMajor[ i ];
  }
  inline       Real& operator () (int i) {
    return symRowMajor[ i ];
  }
};

template<typename Real>
class _SymmetricColMajor : public DataLayout<Real> {
private:
  // Square n-by-n symmetric matrix with (n/2)(n+1) elements 
  // only the upper triangular plus the diagonal are stored.
  int n;
  std::vector<Real> symColMajor; 
  
  // Note: Dont change the order on the indexing integer division 
  // operations (bit shift right) or it might cause underflow.
  // The bit shift right must be the last operation.
  // The indexing for the col-major is a little bit trickier. 
  // math.stackexchange.com/a/2134297
public:
  inline void resize(int numRows_, int numCols_) {
    n = numRows_;
    symColMajor.resize((n * (n + 1)) >> 1);
  }
  inline int rows() const { return n; }
  inline int cols() const { return n; }
  
  inline const Real& operator () (int i, int j) const {
    const int u = j > i ? i : j; // make sure its on upper tri
    const int v = j > i ? j : i;     
    const int index = (((n << 1) * u) - (u * u) -u + (v << 1)) >> 1;
    return symColMajor[index];
  }
  inline       Real& operator () (int i, int j){
    const int u = j > i ? i : j; // make sure its on upper tri
    const int v = j > i ? j : i;     
    const int index = (((n << 1) * u) - (u * u) -u + (v << 1)) >> 1;
    return symColMajor[index];
  }

  // hmm... 
  inline const Real& operator () (int i) const {
    return symColMajor[ i ];
  }
  inline       Real& operator () (int i) {
    return symColMajor[ i ];
  }
};

} // ns Sparse

#endif
