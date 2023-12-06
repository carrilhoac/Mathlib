
#ifndef _Math_EigenSymmetric_h
#define _Math_EigenSymmetric_h 1

template<typename MatrixType>
void Eigen<MatrixType>::tred2() 
{ 
  int q = n-1;

  for (int j = 0; j < n; j++ )
    d[j] = V(q,j);

  for (int i = q; i > 0; i--) {
    DataType scale(0);
    DataType h(0);

    // Scale to avoid under/overflow.
    for (int k = 0; k < i; k++) 
      scale += std::fabs(d[k]);

    if (scale == DataType(0)) {
      e[i] = d[i-1];

      for (int j = 0; j < i; j++) {
        d[j]   = V(i-1,j);
        V(i,j) = DataType(0);
        V(j,i) = DataType(0);
      }
    }
    else {
      // Generate Householder vector.
      for (int k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }

      DataType f = d[i-1];
      DataType g = std::sqrt(h);

      if (f > DataType(0))
        g = -g;

      e[i]   = scale * g;
      h     -= f * g;
      d[i-1] = f - g;

      for (int j = 0; j < i; j++) 
        e[j] = DataType(0);
      
      // Apply similarity transformation to remaining columns.
      for (int j = 0; j < i; j++) {
        f      = d[j];
        V(j,i) = f;
        g      = e[j] + V(j,j) * f;
        for (int k = j+1; k <= i-1; k++) {
          g    += V(k,j) * d[k];
          e[k] += V(k,j) * f;
        }
        e[j] = g;
      }

      f = DataType(0);
      for (int j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      DataType hh = f / (h + h);
      for (int j = 0; j < i; j++) {
        e[j] -= hh * d[j];
      }
      for (int j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (int k = j; k <= i-1; k++) {
            V(k,j) -= (f * e[k] + g * d[k]);
        }
        d[j] = V(i-1,j);
        V(i,j) = DataType(0);
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.
  for (int i = 0; i < q; i++) {
    V(q,i)  = V(i,i);
    V(i,i)  = DataType(1);
    DataType h  = d[i+1];

    if (h != DataType(0)) {
      for (int k = 0; k <= i; k++) 
        d[k] = V(k,i+1) / h;

      for (int j = 0; j <= i; j++) {
        DataType g = DataType(0);

        for (int k = 0; k <= i; k++) 
          g += V(k,i+1) * V(k,j);
        
        for (int k = 0; k <= i; k++) 
          V(k,j) -= g * d[k];
      }
    }
    for (int k = 0; k <= i; k++) 
      V(k,i+1) = DataType(0);
  }
  for (int j = 0; j < n; j++) {
    d[j]   = V(q,j);
    V(q,j) = DataType(0);
  }
  V(q,q) = DataType(1);
  e[0]   = DataType(0);
}


template<typename MatrixType>
void Eigen<MatrixType>::tql2() 
{ 
  for (int i = 1; i < n; i++) 
    e[i-1] = e[i];
  e[n-1] = DataType(0);

  DataType f(0);
  DataType tst1(0);
  //DataType tol = std::pow(2.0,-52.0);
  DataType tol = eps<DataType>();

  for (int l = 0; l < n; l++) {

    // Find small subdiagonal element

    tst1 = std::max(tst1, std::fabs(d[l]) + std::fabs(e[l]));
    int m = l;

    while (m < n) {
      if (std::fabs(e[m]) <= tol*tst1) 
        break;
      m++;
    }

    // If m == l, d[l] is an eigenvalue, otherwise, iterate.

    if (m > l) {
      int iter = 0;
      
      do {
        iter = iter + 1;  // (Could check iteration count here.)

        // Compute implicit shift

        DataType g = d[l];
        DataType p = (d[l+1] - g) / (2.0 * e[l]);
        DataType r = std::hypot(p, DataType(1));

        if (p < 0) 
          r = -r;

        d[l  ] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);

        DataType dl1 = d[l+1];
        DataType h = g - d[l];

        for (int i = l+2; i < n; i++) 
          d[i] -= h;
    
        f = f + h;

        // Implicit QL transformation.

        p = d[m];
        DataType c = DataType(1);
        DataType c2 = c;
        DataType c3 = c;
        DataType el1 = e[l+1];
        DataType s  = DataType(0);
        DataType s2 = DataType(0);

        for (int i = m-1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = std::hypot(p,e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          // Accumulate transformation.

          for (int k = 0; k < n; k++) {
            h = V(k,i+1);
            V(k,i+1) = s * V(k,i) + c * h;
            V(k,i  ) = c * V(k,i) - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        // Check for convergence.

      } while (std::fabs(e[l]) > tol*tst1);
    }
    d[l] = d[l] + f;
    e[l] = DataType(0);
  }
  
  // Sort eigenvalues and corresponding vectors.

  for (int i = 0; i < n-1; i++) {
    int k = i;
    DataType p = d[i];
    for (int j = i+1; j < n; j++) {
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (int j = 0; j < n; j++) {
        p = V(j,i);
        V(j,i) = V(j,k);
        V(j,k) = p;
      }
    }
  }
}

#endif 


