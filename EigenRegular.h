
#ifndef _Math_EigenRegular_h
#define _Math_EigenRegular_h 1

namespace internal 
{
  template<typename Real>
  void complex_div(
    Real xr, Real xi, Real yr, Real yi,  // input
    Real& dr, Real& di)                  // output
  {
    // Complex scalar division
    if (std::fabs(yr) > std::fabs(yi)) {
      Real r = yi / yr;
      Real d = yr + r * yi;
      dr = (xr + r * xi) / d;
      di = (xi - r * xr) / d;
    }
    else {
      Real r = yr / yi;
      Real d = yi + r * yr;
      dr = (r * xr + xi) / d;
      di = (r * xi - xr) / d;
    }
  }
}

template<typename MatrixType>
void Eigen<MatrixType>::orthes(MatrixType& H, VectorType& ort) 
{ 
//  orthes and ortran Fortran subroutines in EISPACK.

  int low = 0;
  int high = n-1;

  for (int m = low+1; m <= high-1; m++) {

    // Scale column.  

    DataType scale(0);
    for (int i = m; i <= high; i++) 
      scale += std::fabs(H(i,m-1));
    
    if (scale != DataType(0)) {

      // Compute Householder transformation.

      DataType h(0);
      for (int i = high; i >= m; i--) {
        ort[i] = H(i,m-1) / scale;
        h += ort[i] * ort[i];
      }

      DataType g = std::sqrt(h);
      if (ort[m] > DataType(0))
        g = -g;
      
      h     -= ort[m] * g;
      ort[m] = ort[m] - g;

      // Apply Householder similarity transformation        
      for (int j = m; j < n; j++) {
        DataType f(0);
        for (int i = high; i >= m; i--) 
          f += ort[i] * H(i,j);
        
        f = f / h;
        for (int i = m; i <= high; i++) 
          H(i,j) -= f * ort[i];
      }

      for (int i = 0; i <= high; i++) {
        DataType f(0);
        for (int j = high; j >= m; j--) 
          f += ort[j] * H(i,j);
        
        f = f / h;
        for (int j = m; j <= high; j++) 
          H(i,j) -= f * ort[j];
      }

      ort[m]   = scale * ort[m];
      H(m,m-1) = scale * g;
    }
  }

  // Accumulate transformations (Algol's ortran).
  V.setIdentity();

  for (int m = high-1; m >= low+1; m--) {
    if ( H(m,m-1) != DataType(0)) {
      for (int i = m+1; i <= high; i++) 
        ort[i] = H(i,m-1);
      
      for (int j = m; j <= high; j++) {

        DataType g = DataType(0.0);
        for (int i = m; i <= high; i++) 
          g += ort[i] * V(i,j);
        
        // Double division avoids possible underflow
        g = (g / ort[m]) / H(m,m-1);
        for (int i = m; i <= high; i++) {
          V(i,j) += g * ort[i];
        }
      }
    }
  }
}

template<typename MatrixType>
void Eigen<MatrixType>::hqr2(MatrixType& H, VectorType& ort) 
{
  int nn = this->n;
  int n = nn-1;
  int low = 0;
  int high = nn-1;

  //Real tol = std::pow(2.0,-52.0);
  DataType tol = eps<DataType>();
  DataType exshift(0);
  DataType p(0),q(0),r(0),s(0),z(0),t,w,x,y;

  // Store roots isolated by balanc and compute matrix norm

  DataType norm(0);
  for (int i = 0; i < nn; i++) {
    if ((i < low) || (i > high)) 
    {
      d[i] = H(i,i);
      e[i] = DataType(0);
    }

    for (int j = std::max(i-1, 0); j < nn; j++) 
      norm += std::fabs(H(i,j));      
  }

  // Outer loop over eigenvalue index

  int iter = 0;
  while (n >= low) {

    // Look for single small sub-diagonal element

    int l = n;
    while (l > low) {
      s = std::fabs(H(l-1,l-1)) + std::fabs(H(l,l));

      if (s == DataType(0)) 
          s = norm;
      
      if (std::fabs(H(l,l-1)) < tol * s) 
          break;
      
      l--;
    }
  
    // Check for convergence
    // One root found

    if (l == n) {
      H(n,n) += exshift;
      d[n] = H(n,n);
      e[n] = DataType(0);
      n--;
      iter = 0;

    // Two roots found

    } else if (l == n-1) {
      w =  H(n  ,n-1) * H(n-1,n);
      p = (H(n-1,n-1) - H(n  ,n)) / DataType(2);
      q = p * p + w;
      z = std::sqrt(std::fabs(q));

      H(n  ,n  ) += exshift;
      H(n-1,n-1) += exshift;
      x = H(n,n);

      // Real pair

      if (q >= 0) {
        if (p >= 0) 
          z = p + z;
        else 
          z = p - z;
        
        d[n-1] = x + z;
        d[n  ] = d[n-1];

        if (z != DataType(0)) 
          d[n] = x - w / z;
        
        e[n-1] = DataType(0);
        e[n  ] = DataType(0);

        x = H(n,n-1);
        s = std::fabs(x) + std::fabs(z);
        p = x / s;
        q = z / s;
        r = std::sqrt(p * p + q * q);
        p = p / r;
        q = q / r;

        // Row modification

        for (int j = n-1; j < nn; j++) {
          z = H(n-1,j);
          H(n-1,j) = q * z + p * H(n,j);
          H(n  ,j) = q * H(n,j) - p * z;
        }

        // Column modification

        for (int i = 0; i <= n; i++) {
          z = H(i,n-1);
          H(i,n-1) = q * z + p * H(i,n);
          H(i,n  ) = q * H(i,n) - p * z;
        }

        // Accumulate transformations

        for (int i = low; i <= high; i++) {
          z = V(i,n-1);
          V(i,n-1) = q * z + p * V(i,n);
          V(i,n  ) = q * V(i,n) - p * z;
        }

      // Complex pair

      } else {
        d[n-1] = x + p;
        d[n  ] = x + p;
        e[n-1] =  z;
        e[n  ] = -z;
      }
      n -= 2;
      iter = 0;

    // No convergence yet

    } else {

      // Form shift

      x = H(n,n);
      y = DataType(0);
      w = DataType(0);

      if (l < n) {
        y = H(n-1,n-1);
        w = H(n  ,n-1) * H(n-1,n  );
      }

      // Wilkinson's original ad hoc shift

      if (iter == 10) {
        exshift += x;
        for (int i = low; i <= n; i++) 
          H(i,i) -= x;
        
        s = std::fabs(H(n,n-1)) + std::fabs(H(n-1,n-2));
        x = y = DataType(0.75) * s;
        w = DataType(-0.4375) * s * s;
      }

      // MATLAB's new ad hoc shift

      if (iter == 30) {
        s = (y - x) / DataType(2);
        s = s * s + w;

        if (s > 0) {
          s = std::sqrt(s);

          if (y < x) 
            s = -s;
          
          s = x - w / ((y - x) / DataType(2) + s);

          for (int i = low; i <= n; i++) 
            H(i,i) -= s;
          
          exshift += s;
          x = y = w = DataType(0.964);
        }
      }

      iter++;   // (Could check iteration count here.)

      // Look for two consecutive small sub-diagonal elements

      int m = n-2;
      while (m >= l) {
        z = H(m,m);
        r = x - z;
        s = y - z;
        p = (r * s - w) / H(m+1,m) + H(m,m+1);

        q = H(m+1,m+1) - z - r - s;
        r = H(m+2,m+1);
        s = std::fabs(p) + std::fabs(q) + std::fabs(r);

        p = p / s;
        q = q / s;
        r = r / s;

        if (m == l) 
          break;
        
        if ((std::fabs(H(m,m-1)) * (std::fabs(q) + std::fabs(r))) <
          tol * (std::fabs(p) * (std::fabs(H(m-1,m-1)) + std::fabs(z) +
          std::fabs(H(m+1,m+1))))) 
        {
          break;
        }
        m--;
      }

      for (int i = m+2; i <= n; i++) {
        H(i,i-2) = DataType(0);

        if (i > m+2) 
          H(i,i-3) = DataType(0);
      }

      // Double QR step involving rows l:n and columns m:n

      for (int k = m; k <= n-1; k++) {
        int notlast = k != n-1;

        if (k != m) {
          p = H(k  ,k-1);
          q = H(k+1,k-1);
          r = notlast ? H(k+2,k-1) : DataType(0);
          x = std::fabs(p) + std::fabs(q) + std::fabs(r);

          if (x != DataType(0)) {
            p = p / x;
            q = q / x;
            r = r / x;
          }
        }

        if (x == DataType(0)) 
          break;

        s = std::sqrt(p*p + q*q + r*r);

        if (p < 0) 
          s = -s;
        
        if (s != 0) {
          if (k != m) 
            H(k,k-1) = -s * x;
          else if (l != m) 
            H(k,k-1) = -H(k,k-1);
            
          p = p + s;
          x = p / s;
          y = q / s;
          z = r / s;
          q = q / p;
          r = r / p;

          // Row modification

          for (int j = k; j < nn; j++) {
            p = H(k,j) + q * H(k+1,j);

            if (notlast) {
              p = p + r * H(k+2,j);
              H(k+2,j)  = H(k+2,j) - p * z;
            }

            H(k  ,j) = H(k  ,j) - p * x;
            H(k+1,j) = H(k+1,j) - p * y;
          }

          // Column modification

          for (int i = 0; i <= std::min(n,k+3); i++) {
            p = x * H(i,k) + y * H(i,k+1);

            if (notlast) {
              p = p + z * H(i,k+2);
              H(i,k+2)  = H(i,k+2) - p * r;
            }

            H(i,k  ) = H(i,k  ) - p;
            H(i,k+1) = H(i,k+1) - p * q;
          }

          // Accumulate transformations

          for (int i = low; i <= high; i++) {
            p = x * V(i,k) + y * V(i,k+1);

            if (notlast) {
              p = p + z * V(i,k+2);
              V(i,k+2)  = V(i,k+2) - p * r;
            }

            V(i,k  ) = V(i,k  ) - p;
            V(i,k+1) = V(i,k+1) - p * q;
          }
        }
      }
    }
  }
  
  // Backsubstitute to find vectors of upper triangular form

  if (norm == DataType(0)) 
    return;

  for (n = nn-1; n >= 0; n--) {
    p = d[n];
    q = e[n];

    // Real vector

    if (q == 0) {
      int l = n;
      H(n,n) = DataType(1);

      for (int i = n-1; i >= 0; i--) {
        w = H(i,i) - p;
        r = DataType(0);
        
        for (int j = l; j <= n; j++) 
          r = r + H(i,j) * H(j,n);
        
        if (e[i] < DataType(0)) {
          z = w;
          s = r;
        } 
        else {
          l = i;

          if (e[i] == DataType(0)) {
            if (w != DataType(0)) 
              H(i,n) = -r / w;
            else 
              H(i,n) = -r / (tol * norm);
            
          // Solve real equations

          } else {
            x = H(i,i+1);
            y = H(i+1,i);
            q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
            t = (x * s - z * r) / q;
            H(i,n) = t;
            if (std::fabs(x) > std::fabs(z)) 
              H(i+1,n) = (-r - w * t) / x;
            else 
              H(i+1,n) = (-s - y * t) / z;
          }

          // Overflow control

          t = std::fabs(H(i,n));
          if ((tol * t) * t > 1) {
            for (int j = i; j <= n; j++) 
              H(j,n) = H(j,n) / t;
          }
        }
      }

    // Complex vector
    } 
    else if (q < 0) {
      int l = n-1;

      // Last vector component imaginary so matrix is triangular

      if (std::fabs(H(n,n-1)) > std::fabs(H(n-1,n))) {
        H(n-1,n-1) = q / H(n,n-1);
        H(n-1,n  ) = -(H(n,n) - p) / H(n,n-1);
      } 
      else {
        internal::complex_div(
          DataType(0), -H(n-1,n), H(n-1,n-1)-p, q, 
          H(n-1,n-1), H(n-1,n));
      }

      H(n,n-1) = DataType(0);
      H(n,n  ) = DataType(1);

      for (int i = n-2; i >= 0; i--) {
        DataType ra,sa,vr,vi;
        ra = DataType(0);
        sa = DataType(0);

        for (int j = l; j <= n; j++) {
          ra = ra + H(i,j) * H(j,n-1);
          sa = sa + H(i,j) * H(j,n  );
        }

        w = H(i,i) - p;

        if (e[i] < DataType(0)) {
          z = w;
          r = ra;
          s = sa;
        } 
        else {
          l = i;

          if (e[i] == 0) {
            internal::complex_div(-ra,-sa,w,q, H(i,n-1), H(i,n));
          } 
          else {

            // Solve complex equations

            x = H(i,i+1);
            y = H(i+1,i);
            vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
            vi = (d[i] - p) * DataType(2) * q;

            if ((vr == DataType(0)) && (vi == DataType(0))) {
              vr = tol * norm * (std::fabs(w) + std::fabs(q) +
              std::fabs(x) + std::fabs(y) + std::fabs(z));
            }

            internal::complex_div(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi, H(i,n-1), H(i,n));

            if (std::fabs(x) > (std::fabs(z) + std::fabs(q))) {
              H(i+1,n-1) = (-ra - w * H(i,n-1) + q * H(i,n  )) / x;
              H(i+1,n  ) = (-sa - w * H(i,n  ) - q * H(i,n-1)) / x;
            } 
            else {
              internal::complex_div(-r-y*H(i,n-1),-s-y*H(i,n),z,q, H(i+1,n-1),H(i+1,n));
            }
          }

          // Overflow control

          t = std::max(std::fabs(H(i,n-1)), std::fabs(H(i,n)));

          if ((tol * t) * t > 1) {
            for (int j = i; j <= n; j++) {
              H(j,n-1) = H(j,n-1) / t;
              H(j,n  ) = H(j,n  ) / t;
            }
          }
        }
      }
    }
  }

  // Vectors of isolated roots

  for (int i = 0; i < nn; i++) {
    if (i < low || i > high) 
      for (int j = i; j < nn; j++) 
        V(i,j) = H(i,j);
  }

  // Back transformation to get eigenvectors of original matrix

  for (int j = nn-1; j >= low; j--) {
    for (int i = low; i <= high; i++) {
      z = DataType(0);

      for (int k = low; k <= std::min(j,high); k++) 
        z = z + V(i,k) * H(k,j);
      
      V(i,j) = z;
    }
  }
}

#endif 

