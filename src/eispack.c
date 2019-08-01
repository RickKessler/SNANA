# include <math.h>
# include <stdbool.h>
# include <stdlib.h>
# include <stdio.h>
# include <time.h>

# include "eispack.h"

/******************************************************************************/

int bakvec ( int n, double t[], double e[], int m, double z[] )

/******************************************************************************/
/*
  Purpose:

    BAKVEC determines eigenvectors by reversing the FIGI transformation.

  Discussion:

    BAKVEC forms the eigenvectors of a nonsymmetric tridiagonal
    matrix by back transforming those of the corresponding symmetric
    matrix determined by FIGI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double T[N*3], contains the nonsymmetric matrix.  Its
    subdiagonal is stored in the positions 2:N of the first column,
    its diagonal in positions 1:N of the second column,
    and its superdiagonal in positions 1:N-1 of the third column.
    T(1,1) and T(N,3) are arbitrary.

    Input/output, double E[N].  On input, E(2:N) contains the
    subdiagonal elements of the symmetric matrix.  E(1) is arbitrary.
    On output, the contents of E have been destroyed.

    Input, int M, the number of eigenvectors to be back
    transformed.

    Input/output, double Z[N*M], contains the eigenvectors.
    On output, they have been transformed as requested.

    Output, int BAKVEC, an error flag.
    0, for normal return,
    2*N+I, if E(I) is zero with T(I,1) or T(I-1,3) non-zero.
    In this case, the symmetric matrix is not similar
    to the original matrix, and the eigenvectors
    cannot be found by this program.
*/
{
  int i;
  int ierr;
  int j;

  ierr = 0;

  if ( m == 0 )
  {
    return ierr;
  }

  e[0] = 1.0;
  if ( n == 1 )
  {
    return ierr;
  }

  for ( i = 1; i < n; i++ )
  {
    if ( e[i] == 0.0 )
    {
      if ( t[i+0*3] != 0.0 || t[i-1+2*3] != 0.0 )
      {
        ierr = 2 * n + ( i + 1 );
        return ierr;
      }
      e[i] = 1.0;
    }
    else
    {
      e[i] = e[i-1] * e[i] / t[i-1+2*3];
    }
  }

  for ( j = 0; j < m; j++ )
  {
    for ( i = 1; i < n; i++ )
    {
      z[i+j*n] = z[i+j*n] * e[i];
    }
  }

  return ierr;
}
/******************************************************************************/

void balanc ( int n, double a[], int *low, int *igh, double scale[] )

/******************************************************************************/
/*
  Purpose:

    BALANC balances a real matrix before eigenvalue calculations.

  Discussion:

    BALANC balances a real matrix and isolates eigenvalues.

    Suppose that the principal submatrix in rows LOW through IGH
    has been balanced, that P(J) denotes the index interchanged
    with J during the permutation step, and that the elements
    of the diagonal matrix used are denoted by D(I,J).  Then

      SCALE(J) = P(J),    J = 1,...,LOW-1,
               = D(J,J),  J = LOW,...,IGH,
               = P(J)     J = IGH+1,...,N.

    The order in which the interchanges are made is N to IGH+1,
    then 1 to LOW-1.

    Note that 1 is returned for LOW if IGH is zero formally.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double A(N,N), the N by N matrix.  On output,
    the matrix has been balanced.

    Output, int LOW, IGH, indicate that A(I,J) is equal to
    zero if
    (1) I is greater than J and
    (2) J=1,...,LOW-1 or I=IGH+1,...,N.

    Output, double SCALE(N), contains information determining the
    permutations and scaling factors used.
*/
{
  double b2;
  double c;
  bool done;
  double f;
  double g;
  int i;
  int j;
  int k;
  int l;
  int m;
  bool noconv;
  double r;
  double radix;
  double s;
  bool swap;
  double t;

  radix = 16.0;
  b2 = radix * radix;
  j = -1;
  m = -1;
  k = 0;
  l = n - 1;
/*
  Search for rows isolating an eigenvalue and push them down.
*/
  done = false;

  while ( ! done )
  {
    for ( j = l; 0 <= j; j-- )
    {
      swap = true;
      for ( i = 0; i <= l; i++ )
      {
        if ( i != j )
        {
          if ( a[j+i*n] != 0.0 )
          {
            swap = false;
            break;
          }
        }
      }

      if ( swap )
      {
        m = l;
        scale[m] = ( double ) j;

        if ( j != m )
        {
          for ( i = 0; i <= l; i++ )
          {
            t        = a[i+j*n];
            a[i+j*n] = a[i+m*n];
            a[i+m*n] = t;
          }

          for ( i = k; i < n; i++ )
          {
            t        = a[j+i*n];
            a[j+i*n] = a[m+i*n];
            a[m+i*n] = t;
          }
        }

        if ( l == 0 )
        {
          *low = k;
          *igh = l;
          return;
        }

        l = l - 1;
        if ( l < 0 )
        {
          done = true;
        }
        break;
      }
      else if ( j == 0 )
      {
        done = true;
        break;
      }
    }
  }
/*
  Search for columns isolating an eigenvalue and push them left.
*/
  done = false;

  while ( ! done )
  {
    for ( j = k; j <= l; j++ )
    {
      swap = true;

      for ( i = k; i <= l; i++ )
      {
        if ( i != j )
        {
          if ( a[i+j*n] != 0.0 )
          {
            swap = false;
            break;
          }
        }
      }

      if ( swap )
      {
        m = k;
        scale[m] = ( double ) j;

        if ( j != m )
        {
          for ( i = 0; i <= l; i++ )
          {
            t        = a[i+j*n];
            a[i+j*n] = a[i+m*n];
            a[i+m*n] = t;
          }
          for ( i = k; i < n; i++ )
          {
            t        = a[j+i*n];
            a[j+i*n] = a[m+i*n];
            a[m+i*n] = t;
          }
        }

        k = k + 1;
        if ( l < k )
        {
          done = true;
        }
        break;
      }
      else
      {
        if ( j == l )
        {
          done = true;
          break;
        }
      }
    }
  }
/*
  Balance the submatrix in rows K to L.
*/
  for ( i = k; i <= l; i++ )
  {
    scale[i] = 1.0;
  }
/*
  Iterative loop for norm reduction.
*/
  noconv = true;

  while ( noconv )
  {
    noconv = false;

    for ( i = k; i <= l; i++ )
    {
      c = 0.0;
      r = 0.0;
      for ( j = k; j <= l; j++ )
      {
        if ( j != i )
        {
          c = c + fabs ( a[j+i*n] );
          r = r + fabs ( a[i+j*n] );
        }
      }
/*
  Guard against zero C or R due to underflow.
*/
      if ( c != 0.0 && r != 0.0 )
      {
        g = r / radix;
        f = 1.0;
        s = c + r;

        while ( c < g )
        {
          f = f * radix;
          c = c * b2;
        }

        g = r * radix;
        while ( g <= c )
        {
          f = f / radix;
          c = c / b2;
        }
/*
  Balance.
*/
        if ( ( c + r ) / f < 0.95 * s )
        {
          g = 1.0 / f;
          scale[i] = scale[i] * f;
          noconv = true;

          for ( j = k; j < n; j++ )
          {
            a[i+j*n] = a[i+j*n] * g;
          }
          for ( j = 0; j <= l; j++ )
          {
            a[j+i*n] = a[j+i*n] * f;
          }
        }
      }
    }
  }

  *low = k;
  *igh = l;

  return;
}
/******************************************************************************/

void balbak ( int n, int low, int igh, double scale[], int m, double z[] )

/******************************************************************************/
/*
  Purpose:

    BALBAK determines eigenvectors by undoing the BALANC transformation.

  Discussion:

    BALBAK forms the eigenvectors of a real general matrix by
    back transforming those of the corresponding balanced matrix
    determined by BALANC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 July 2013

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Parlett and Reinsch,
    Numerische Mathematik,
    Volume 13, pages 293-304, 1969.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, column indices determined by BALANC.

    Input, double SCALE[N], contains information determining
    the permutations and scaling factors used by BALANC.

    Input, int M, the number of columns of Z to be
    back-transformed.

    Input/output, double Z[N*M], contains the real and imaginary 
    parts of the eigenvectors, which, on return, have been back-transformed.
*/
{
  int i;
  int ii;
  int j;
  int k;
  double t;

  if ( m <= 0 )
  {
    return;
  }

  if ( igh != low )
  {
    for ( i = low; i <= igh; i++ )
    {
      for ( j = 0; j < m; j++ )
      {
        z[i+j*n] = z[i+j*n] * scale[i];
      }
    }
  }

  for ( ii = 0; ii < n; ii++ )
  {
    i = ii;

    if ( i < low || igh < i )
    {
      if ( i < low )
      {
        i = low - ii;
      }

      k = ( int ) ( scale[i] );

      if ( k != i )
      {
        for ( j = 0; j < m; j++ )
        {
          t        = z[i+j*n];
          z[i+j*n] = z[k+j*n];
          z[k+j*n] = t;
        }
      }
    }
  }

  return;
}
/******************************************************************************/

void bandr ( int n, int mb, double a[], double d[], double e[], double e2[], 
  bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    BANDR reduces a symmetric band matrix to symmetric tridiagonal form.

  Discussion:

    BANDR reduces a real symmetric band matrix
    to a symmetric tridiagonal matrix using and optionally
    accumulating orthogonal similarity transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int MB, is the (half) band width of the matrix,
    defined as the number of adjacent diagonals, including the principal
    diagonal, required to specify the non-zero portion of the
    lower triangle of the matrix.

    Input/output, double A[N*MB].  On input, contains the lower
    triangle of the symmetric band input matrix stored as an N by MB array.
    Its lowest subdiagonal is stored in the last N+1-MB positions of the first
    column, its next subdiagonal in the last N+2-MB positions of the second
    column, further subdiagonals similarly, and finally its principal diagonal
    in the N positions of the last column.  Contents of storages not part of
    the matrix are arbitrary.  On output, A has been destroyed, except for
    its last two columns which contain a copy of the tridiagonal matrix.

    Output, double D[N], the diagonal elements of the tridiagonal
    matrix.

    Output, double E[N], the subdiagonal elements of the tridiagonal
    matrix in E(2:N).  E(1) is set to zero.

    Output, double E2[N], contains the squares of the corresponding
    elements of E.  E2 may coincide with E if the squares are not needed.

    Input, bool MATZ, is true if the transformation matrix is
    to be accumulated, and false otherwise.

    Output, double Z[N*N], the orthogonal transformation matrix
    produced in the reduction if MATZ is true.  Otherwise, Z is
    not referenced.
*/
{
  double b1;
  double b2;
  double c2;
  double dmin;
  double dminrt;
  double f1;
  double f2;
  double g;
  int i;
  int i1;
  int i2;
  int j;
  int j1;
  int j2;
  int jj;
  int k;
  int kr;
  int l;
  int m1;
  int maxl;
  int maxr;
  int mr;
  int r;
  int r1;
  double s2;
  double u;
  int ugl;

  dmin = r8_epsilon ( );
  dminrt = sqrt ( dmin );
/*
  Initialize the diagonal scaling matrix.
*/
  for ( i = 0; i < n; i++ )
  {
    d[i] = 1.0;
  }

  if ( matz )
  {
    r8mat_identity ( n, z );
  }
/*
  Is input matrix diagonal?
*/
  if ( mb == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      d[i] = a[i+(mb-1)*n];
      e[i] = 0.0;
      e2[i] = 0.0;
    }
    return;
  }

  m1 = mb - 1;

  if ( m1 != 1 )
  {
    for ( k = 1; k <= n - 2; k++ )
    {
      maxr = i4_min ( m1, n - k );
      for ( r1 = 2; r1 <= maxr; r1++ )
      {
        r = maxr + 2 - r1;
        kr = k + r;
        mr = mb - r;
        g = a[kr-1+(mr-1)*n];
        a[kr-2+0*n] = a[kr-2+mr*n];
        ugl = k;

        for ( j = kr; j <= n; j = j + m1 )
        {
          j1 = j - 1;
          j2 = j1 - 1;

          if ( g == 0.0 )
          {
            break;
          }

          b1 = a[j1-1+0*n] / g;
          b2 = b1 * d[j1-1] / d[j-1];
          s2 = 1.0 / ( 1.0 + b1 * b2 );

          if ( s2 < 0.5 )
          {
            b1 = g / a[j1-1+0*n];
            b2 = b1 * d[j-1] / d[j1-1];
            c2 = 1.0 - s2;
            d[j1-1] = c2 * d[j1-1];
            d[j-1] = c2 * d[j-1];
            f1 = 2.0 * a[j-1+(m1-1)*n];
            f2 = b1 * a[j1-1+(mb-1)*n];
            a[j-1+(m1-1)*n] = - b2 * ( b1 * a[j-1+(m1-1)*n] - a[j-1+(mb-1)*n] ) 
              - f2 + a[j-1+(m1-1)*n];
            a[j1-1+(mb-1)*n] = b2 * ( b2 * a[j-1+(mb-1)*n] + f1 ) + a[j1-1+(mb-1)*n];
            a[j-1+(mb-1)*n] = b1 * ( f2 - f1 ) + a[j-1+(mb-1)*n];

            for ( l = ugl; l <= j2; l++ )
            {
              i2 = mb - j + l;
              u = a[j1-1+i2*n] + b2 * a[j-1+(i2-1)*n];
              a[j-1+(i2-1)*n] = -b1 * a[j1-1+i2*n] + a[j-1+(i2-1)*n];
              a[j1-1+i2*n] = u;
            }

            ugl = j;
            a[j1-1+0*n] = a[j1+0*n] + b2 * g;

            if ( j != n )
            {
              maxl = i4_min ( m1, n - j1 );

              for ( l = 2; l <= maxl; l++ )
              {
                i1 = j1 + l;
                i2 = mb - l;
                u = a[i1-1+(i2-1)*n] + b2 * a[i1-1+i2*n];
                a[i1-1+i2*n] = -b1 * a[i1-1+(i2-1)*n] + a[i1-1+i2*n];
                a[i1-1+(i2-1)*n] = u;
              }

              i1 = j + m1;

              if ( i1 <= n )
              {
                g = b2 * a[i1-1+0*n];
              }
            }

            if ( matz )
            {
              for ( l = 1; l <= n; l++ )
              {
                u = z[l-1+(j1-1)*n] + b2 * z[l-1+(j-1)*n];
                z[l-1+(j-1)*n] = -b1 * z[l-1+(j1-1)*n] + z[l-1+(j-1)*n];
                z[l-1+(j1-1)*n] = u;
              }
            }
          }
          else
          {
            u = d[j1-1];
            d[j1-1] = s2 * d[j-1];
            d[j-1] = s2 * u;
            f1 = 2.0 * a[j-1+(m1-1)*n];
            f2 = b1 * a[j-1+(mb-1)*n];
            u = b1 * ( f2 - f1 ) + a[j1-1+(mb-1)*n];
            a[j-1+(m1-1)*n] = b2 * ( b1 * a[j-1+(m1-1)*n] - a[j1-1+(mb-1)*n] ) 
              + f2 - a[j-1+(m1-1)*n];
            a[j1-1+(mb-1)*n] = b2 * ( b2 * a[j1-1+(mb-1)*n] + f1 ) 
              + a[j-1+(mb-1)*n];
            a[j-1+(mb-1)*n] = u;

            for ( l = ugl; l <= j2; l++ )
            {
              i2 = mb - j + l;
              u = b2 * a[j1-1+i2*n] + a[j-1+(i2-1)*n];
              a[j-1+(i2-1)*n] = -a[j1-1+i2*n] + b1 * a[j-1+(i2-1)*n];
              a[j1-1+i2*n] = u;
            }

            ugl = j;
            a[j1-1+0*n] = b2 * a[j1-1+0*n] + g;
 
            if ( j != n )
            {
              maxl = i4_min ( m1, n - j1 );

              for ( l = 2; l <= maxl; l++ )
              {
                i1 = j1 + l;
                i2 = mb - l;
                u = b2 * a[i1-1+(i2-1)*n] + a[i1-1+i2*n];
                a[i1-1+i2*n] = -a[i1-1+(i2-1)*n] + b1 * a[i1-1+i2*n];
                a[i1-1+(i2-1)*n] = u;
              }

              i1 = j + m1;

              if ( i1 <= n )
              {
                g = a[i1-1+0*n];
                a[i1-1+0*n] = b1 * a[i1-1+0*n];
              }
            }

            if ( matz )
            {
              for ( l = 1; l <= n; l++ )
              {
                u = b2 * z[l-1+(j1-1)*n] + z[l-1+(j-1)*n];
                z[l-1+(j-1)*n] = -z[l-1+(j1-1)*n] + b1 * z[l-1+(j-1)*n];
                z[l-1+(j1-1)*n] = u;
              }
            }
          }
        }
      }
/*
  Rescale to avoid underflow or overflow.
*/
      if ( ( k % 64 ) == 0 )
      {
        for ( j = k; j <= n; j++ )
        {
          if ( d[j-1] < dmin )
          {
            maxl = i4_max ( 1, mb + 1 - j );

            for ( jj = maxl; jj <= m1; jj++ )
            {
              a[j-1+(jj-1)*n] = dminrt * a[j-1+(jj-1)*n];
            }

            if ( j != n )
            {
              maxl = i4_min ( m1, n - j );

              for ( l = 1; l <= maxl; l++ )
              {
                i1 = j + l;
                i2 = mb - l;
                a[i1-1+(i2-1)*n] = dminrt * a[i1-1+(i2-1)*n];
              }
            }

            if ( matz )
            {
              for ( i = 1; i <= n; i++ )
              {
                z[i-1+(j-1)*n] = dminrt * z[i-1+(j-1)*n];
              }
            }

            a[j-1+(mb-1)*n] = dmin * a[j-1+(mb-1)*n];
            d[j-1] = d[j-1] / dmin;
          }
        }
      }
    }
  }
/*
  Form square root of scaling matrix.
*/
  for ( i = 1; i < n; i++ )
  {
    e[i] = sqrt ( d[i] );
  }
  if ( matz )
  {
    for ( j = 1; j < n; j++ )
    {
      for ( i = 0; i < n; i++ )
      {
        z[i+j*n] = z[i+j*n] * e[j];
      }
    }
  }

  u = 1.0;

  for ( j = 1; j < n; j++ )
  {
    a[j+(m1-1)*n] = u * e[j] * a[j+(m1-1)*n];
    u = e[j];
    e2[j] = a[j+(m1-1)*n] * a[j+(m1-1)*n];
    a[j+(mb-1)*n] = d[j] * a[j+(mb-1)*n];
    d[j] = a[j+(mb-1)*n];
    e[j] = a[j+(m1-1)*n];
  }

  d[0] = a[0+(mb-1)*n];
  e[0] = 0.0;
  e2[0] = 0.0;

  return;
}
/******************************************************************************/

int bandv ( int n, int mbw, double a[], double e21, int m, double w[], 
  double z[] )

/******************************************************************************/
/*
  Purpose:

    BANDV finds eigenvectors from eigenvalues, for a real symmetric band matrix.

  Discussion:

    BANDV finds those eigenvectors of a real symmetric
    band matrix corresponding to specified eigenvalues, using inverse
    iteration.  

    The routine may also be used to solve systems of linear equations with a 
    symmetric or non-symmetric band coefficient matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int MBW, the number of columns of the array A used
    to store the band matrix.  If the matrix is symmetric, MBW is its (half)
    band width, denoted MB and defined as the number of adjacent
    diagonals, including the principal diagonal, required to
    specify the non-zero portion of the lower triangle of the
    matrix.  If the routine is being used to solve systems
    of linear equations and the coefficient matrix is not
    symmetric, it must however have the same number of adjacent
    diagonals above the main diagonal as below, and in this
    case, MBW=2*MB-1.

    Input, double A[N,MBW], the lower triangle of the symmetric
    band input matrix stored as an N by MB array.  Its lowest subdiagonal is
    stored in the last N+1-MB positions of the first column, its next
    subdiagonal in the last N+2-MB positions of the second column, further
    subdiagonals similarly, and finally its principal diagonal in the N
    positions of column MB.  If the routine is being used to solve systems
    of linear equations, and the coefficient matrix is not symmetric, A is
    N by 2*MB-1 instead, with lower triangle as above and with its first
    superdiagonal stored in the first N-1 positions of column MB+1, its
    second superdiagonal in the first N-2 positions of column MB+2, further
    superdiagonals similarly, and finally its highest superdiagonal in
    the first N+1-MB positions of the last column.  Contents of storages
    not part of the matrix are arbitrary.

    Input, double E21, specifies the ordering of the eigenvalues
    and contains 0.0 if the eigenvalues are in ascending order, or 2.0 if
    the eigenvalues are in descending order.  If the routine is being used
    to solve systems of linear equations, E21 should be set to 1.0
    if the coefficient matrix is symmetric and to -1.0 if not.

    Input, int M, the number of specified eigenvalues or the
    number of systems of linear equations.

    Input, double W[M], contains the M eigenvalues in ascending or
    descending order.  If the routine is being used to solve systems of
    linear equations (A-W(1:M)*I) * X(1:M) = B(1:M), where I is the identity
    matrix, W should be set accordingly.

    Input/output, double Z[N,M].  On input, the constant matrix
    columns B(1:M), if the routine is used to solve systems of linear
    equations.  On output, the associated set of orthogonal eigenvectors.
    Any vector which fails to converge is set to zero.  If the
    routine is used to solve systems of linear equations,
    Z contains the solution matrix columns X(1:M).

    Output, int BANDV, error flag.
    0, for normal return,
    -R, if the eigenvector corresponding to the R-th eigenvalue fails to
    converge, or if the R-th system of linear equations is nearly singular.
*/
{
  double eps2 = 0.0 ;
  double eps3 = 0.0 ;
  double eps4 = 0.0 ;
  int group = 0;
  int i;
  int ierr;
  int ii;
  int ij;
  int ij1;
  int its;
  int j;
  int jj;
  int k;
  int kj;
  int kj1;
  int m1;
  int m21;
  int maxj;
  int maxk;
  int mb;
  double norm;
  double order;
  int r;
  double *rv;
  double *rv6;
  double t;
  double u;
  double uk = 0.0 ;
  double v;
  double x0;
  double x1;
  double xu;

  ierr = 0;

  if ( m == 0 )
  {
    return ierr;
  }

  rv = ( double * ) malloc ( n * ( 2 * mbw - 1 ) * sizeof ( double ) );
  rv6 = ( double * ) malloc ( n * sizeof ( double ) );

  x0 = 0.0;

  if ( e21 < 0.0 )
  {
    mb = ( mbw + 1 ) / 2;
  }
  else
  {
    mb = mbw;
  }

  m1 = mb - 1;
  m21 = m1 + mb;
  order = 1.0 - fabs ( e21 );
/*
  Find vectors by inverse iteration.
*/
  for ( r = 1; r <= m; r++ )
  {
    its = 1;
    x1 = w[r-1];
/*
  Compute norm of matrix.
*/
    if ( r == 1 )
    {
      norm = 0.0;

      for ( j = 1; j <= mb; j++ )
      {
        jj = mb + 1 - j;
        kj = jj + m1;
        ij = 1;

        v = 0.0;
        for ( i = mb + 1 - j; i <= n; i++ )
        {
          v = v + fabs ( a[i-1+(j-1)*n] );
          if ( e21 < 0.0 )
          {
            v = v + fabs ( a[ij-1+(kj-1)*n] );
            ij = ij + 1;
          }
        }
        norm = r8_max ( norm, v );
      }

      if ( e21 < 0.0 )
      {
        norm = 0.5 * norm;
      }
/*
  EPS2 is the criterion for grouping,
  EPS3 replaces zero pivots and equal roots are modified by eps3,
  EPS4 is taken very small to avoid overflow.
*/
      if ( norm == 0.0 )
      {
        norm = 1.0;
      }

      eps2 = 0.001 * norm * fabs ( order);
      eps3 = fabs ( norm ) * r8_epsilon ( );
      uk = n;
      uk = sqrt ( uk );
      eps4 = uk * eps3;
      group = 0;
    }
/*
  Look for close or coincident roots.
*/
    else
    {
      if ( eps2 <= fabs ( x1 - x0 ) )
      {
        group = 0;
      }
      else
      {
        group = group + 1;

        if ( order * ( x1 - x0 ) <= 0.0 )
        {
          x1 = x0 + order * eps3;
        }
      }
    }
/*
  Expand matrix, subtract eigenvalue, and initialize vector.
*/
    for ( i = 1; i <= n; i++ )
    {
      ij = i + i4_min ( 0, i - m1 ) * n;
      kj = ij + mb * n;
      ij1 = kj + m1 * n;

      for ( j = 1; j <= m1; j++ )
      {
        if ( ij <= m1 )
        {
          if ( ij <= 0 )
          {
            rv[ij1-1] = 0.0;
            ij1 = ij1 + n;
          }
        }
        else
        {
          rv[ij-1] = a[i-1+(j-1)*n];
        }

        ij = ij + n;
        ii = i + j;

        if ( ii <= n )
        {
          jj = mb - j;

          if ( e21 < 0.0 )
          {
            ii = i;
            jj = mb + j;
          }
          rv[kj-1] = a[ii-1+(jj-1)*n];
          kj = kj + n;
        }
      }

      rv[ij-1] = a[i-1+(mb-1)*n] - x1;
      rv6[i-1] = eps4;
      if ( order == 0.0 )
      {
        rv6[i-1] = z[i-1+(r-1)*n];
      }
    }

    if ( m1 != 0 )
    {
/*
  Elimination with interchanges.
*/
      for ( i = 1; i <= n; i++ )
      {
        ii = i + 1;
        maxk = i4_min ( i + m1 - 1, n );
        maxj = i4_min ( n - i, m21 - 2 ) * n;

        for ( k = i; k <= maxk; k++ )
        {
          kj1 = k;
          j = kj1 + n;
          jj = j + maxj;

          for ( kj = j; kj <= jj; kj = kj + n )
          {
            rv[kj1-1] = rv[kj-1];
            kj1 = kj;
          }
          rv[kj1-1] = 0.0;
        }

        if ( i < n )
        {
          u = 0.0;
          maxk = i4_min ( i + m1, n );
          maxj = i4_min ( n - ii, m21 - 2 ) * n;

          for ( j = i; j <= maxk; j++ )
          {
            if ( fabs ( u ) <= fabs ( rv[j-1] ) )
            {
              u = rv[j-1];
              k = j;
            }
          }

          j = i + n;
          jj = j + maxj;

          if ( k != i )
          {
            kj = k;

            for ( ij = i; ij <= jj; ij = ij + n )
            {
              t        = rv[ij-1];
              rv[ij-1] = rv[kj-1];
              rv[kj-1] = t;
              kj = kj + n;
            }

            if ( order == 0.0 )
            {
              t        = rv6[i-1];
              rv6[i-1] = rv6[k-1];
              rv6[k-1] = t;
            }
          }

          if ( u != 0.0 )
          {
            for ( k = ii; k <= maxk; k++ )
            {
              v = rv[k-1] / u;
              kj = k;

              for ( ij = j; ij <= jj; ij = ij + n )
              {
                kj = kj + n;
                rv[kj-1] = rv[kj-1] - v * rv[ij-1];
              }

              if ( order == 0.0 )
              {
                rv6[k-1] = rv6[k-1] - v * rv6[i-1];
              }
            }
          }
        }
      }
    }
/*
  Back substitution.
*/
    while ( true )
    {
      for ( i = n; 1 <= i; i-- )
      {
        maxj = i4_min ( n + 1 - i, m21 );

        if ( maxj != 1 )
        {
          ij1 = i;
          j = ij1 + n;
          jj = j + ( maxj - 2 ) * n;

          for ( ij = j; ij <= jj; ij = ij + n )
          {
            ij1 = ij1 + 1;
            rv6[i-1] = rv6[i-1] - rv[ij-1] * rv6[ij1-1];
          }
        }

        v = rv[i-1];
/*
  Error: nearly singular linear system.
*/
        if ( fabs ( v ) < eps3 )
        {
          if ( order == 0.0 )
          {
            ierr = - r;
          }
          v = fabs ( eps3 ) * r8_sign ( v );
        }
        rv6[i-1] = rv6[i-1] / v;
      }

      xu = 1.0;

      if ( order == 0.0 )
      {
        for ( j = 1; j <= n; j++ )
        {
          z[j-1+(r-1)*n] = rv6[j-1] * xu;
        }
        x0 = x1;
        break;
      }
/*
  Orthogonalize with respect to previous members of group.
*/
      for ( j = r - group; j <= r - 1; j++ )
      {
        xu = 0.0;
        for ( i = 1; i <= n; i++ )
        {
          xu = xu + rv6[i-1] * z[i-1+(j-1)*n];
        }
        for ( i = 1; i <= n; i++ )
        {
          rv6[i-1] = rv6[i-1] - xu * z[i-1+(j-1)*n];
        }
      }

      norm = 0.0;
      for ( i = 1; i <= n; i++ )
      {
        norm = norm + fabs ( rv6[i-1] );
      }
/*
  Choose a new starting vector.
*/
      if ( 0.1 <= norm )
      {
        u = 0.0;
        for ( i = 1; i <= n; i++ )
        {
          u = pythag ( u, rv6[i-1] );
        }
        xu = 1.0 / u;
        for ( i = 1; i <= n; i++ )
        {
          z[i-1+(r-1)*n] = rv6[i-1] * xu;
        }
        x0 = x1;
        break;
      }
      else if ( n <= its )
      {
        ierr = - r;
        xu = 0.0;
        for ( i = 1; i <= n; i++ )
        {
          z[i-1+(r-1)*n] = rv6[i-1] * xu;
        }
        x0 = x1;
        break;
      }
      else
      {
        its = its + 1;
        xu = eps4 / ( uk + 1.0 );
        rv6[0] = eps4;
        for ( i = 2; i <= n; i++ )
        {
          rv6[i-1] = xu;
        }
        rv6[its-1] = rv6[its-1] - eps4 * uk;
      }
    }
  }

  free ( rv );
  free ( rv6 );

  return ierr;
}
/******************************************************************************/

int bisect ( int n, double *eps1, double d[], double e[], double e2[], 
  double t1, double t2, int mm, int *m, double w[], int ind[] )

/******************************************************************************/
/*
  Purpose:

    BISECT computes some eigenvalues of a real symmetric tridiagonal matrix.

  Discussion:

    BISECT finds those eigenvalues of a real symmetric tridiagonal matrix 
    which lie in a specified interval, using bisection.

    In the original code, the lower and upper seach bounds were
    copied, then modified, and then restored at the end.  But they
    should really be input only

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 March 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double *EPS1, is an absolute error tolerance for
    the computed eigenvalues.  If the input EPS1 is non-positive, it is reset
    for each submatrix to a default value, namely, minus the product of the
    relative machine precision and the 1-norm of the submatrix.

    Input, double D[N], the diagonal elements of the input matrix.

    Input, double E[N], contains in E(2:N) the subdiagonal elements
    of the matrix.  E[0] is arbitrary.

    Input/output, double E2[N].  On input, the squares of the
    corresponding elements of E.  E2(1) is arbitrary.  On output, elements of 
    E2, corresponding to elements of E regarded as negligible, have been
    replaced by zero, causing the matrix to split into a direct sum of
    submatrices.  E2[0] is also set to zero.

    Input, double T1, T2, define the interval to be searched for
    eigenvalues.  If T1 is not less than T2, no eigenvalues will be found.

    Input, int MM, an upper bound for the number of
    eigenvalues in the interval.  Warning: if more than MM eigenvalues are
    determined to lie in the interval, an error return is made with no
    eigenvalues found.

    Output, int M, the number of eigenvalues determined to lie
    in (LB,UB).

    Output, double W[M], the eigenvalues in ascending order.

    Output, int IND[MM], contains in its first M positions
    the submatrix indices associated with the corresponding eigenvalues in W:
    0 for eigenvalues belonging to the first submatrix from the top, 1 for
    those belonging to the second submatrix, and so on.

    Output, int BISECT, error flag.
    0, for normal return,
    3*N+1, if M exceeds MM.
*/
{
  int i;
  int ierr;
  int j;
  int k;
  int l;
  double lb;
  int m1;
  int m2;
  int p;
  int q;
  int r;
  double *rv4;
  double *rv5;
  int s;
  int tag;
  double tst1;
  double tst2;
  double u;
  double ub;
  double v;
  double x0;
  double x1;
  double xu;

  rv4 = ( double * ) malloc ( n * sizeof ( double ) );
  rv5 = ( double * ) malloc ( n * sizeof ( double ) );

  ierr = 0;
  s = 0;
  tag = 0;
  lb = t1;
  ub = t2;
/*
  Look for small sub-diagonal entries.
*/
  e2[0] = 0.0;

  for ( i = 2; i <= n; i++ )
  {
    tst1 = fabs ( d[i-1] ) + fabs ( d[i-2] );
    tst2 = tst1 + fabs ( e[i-1] );

    if ( tst2 <= tst1 )
    {
      e2[i-1] = 0.0;
    }
  }
/*
  Determine the number of eigenvalues in the interval.
*/
  p = 1;
  q = n;

  x1 = ub;
  s = sturm_sequence ( d, e, e2, n, p - 1, q - 1, x1 );
  *m = s;

  x1 = lb;
  s = sturm_sequence ( d, e, e2, n, p - 1, q - 1, x1 );
  *m = *m - s;

  if ( mm < *m )
  {
    ierr = 3 * n + 1;
    free ( rv4 );
    free ( rv5 );
    return ierr;
  }

  q = 0;
  r = 0;
/*
  Establish and process next submatrix, refining
  interval by the Gerschgorin bounds.
*/
  while ( true )
  {
    if ( r == *m )
    {
      free ( rv4 );
      free ( rv5 );
      return ierr;
    }
    tag = tag + 1;
    p = q + 1;
    xu = d[p-1];
    x0 = d[p-1];
    u = 0.0;

    for ( q = p; q <= n; q++ )
    {
      x1 = u;
      u = 0.0;
      v = 0.0;

      if ( q < n )
      {
        u = fabs ( e[q] );
        v = e2[q];
      }
      xu = r8_min ( d[q-1] - ( x1 + u ), xu );
      x0 = r8_max ( d[q-1] + ( x1 + u ), x0 );

      if ( v == 0.0 )
      {
        break;
      }
    }

    x1 = r8_max ( fabs ( xu ), fabs ( x0 ) ) * r8_epsilon ( );
    if ( *eps1 <= 0.0 )
    {
      *eps1 = - x1;
    }
/*
  Check for an isolated root within interval.
*/
    if ( p == q )
    {
      if ( d[p-1] < t1 || t2 <= d[p-1] )
      {
        if ( q < n )
        {
          continue;
        }
        else
        {
          free ( rv4 );
          free ( rv5 );
          return ierr;
        }
      }
      m1 = p;
      m2 = p;
      rv5[p-1] = d[p-1];
    }
    else
    {
      x1 = x1 * ( q - p + 1 );
      lb = r8_max ( t1, xu - x1 );
      ub = r8_min ( t2, x0 + x1 );
      x1 = lb;
      s = sturm_sequence ( d, e, e2, n, p - 1, q - 1, x1 );
      m1 = s + 1;
      x1 = ub;
      s = sturm_sequence ( d, e, e2, n, p - 1, q - 1, x1 );
      m2 = s;

      if ( m2 < m1 )
      {
        if ( q < n )
        {
          continue;
        }
        else
        {
          free ( rv4 );
          free ( rv5 );
          return ierr;
        }
      }
/*
  Find roots by bisection.
*/
      x0 = ub;
      for ( i = m1; i <= m2; i++ )
      {
        rv5[i-1] = ub;
        rv4[i-1] = lb;
      }
/*
  Loop for the K-th eigenvalue.
*/
      k = m2;

      while ( true )
      {
        xu = lb;
        for ( i = k; m1 <= i; i-- )
        {
          if ( xu < rv4[i-1] )
          {
            xu = rv4[i-1];
            break;
          }
        }
        x0 = r8_min ( x0, rv5[k-1] );
/*
  Next bisection step.
*/
        while ( true )
        {
          x1 = ( xu + x0 ) * 0.5;

          if ( ( x0 - xu ) <= fabs ( *eps1 ) )
          {
            break;
          }
          tst1 = 2.0 * ( fabs ( xu ) + fabs ( x0 ) );
          tst2 = tst1 + ( x0 - xu );
          if ( tst2 == tst1 )
          {
            break;
          }
          s = sturm_sequence ( d, e, e2, n, p - 1, q - 1, x1 );
/*
  Refine intervals.
*/
          if ( k <= s )
          {
            x0 = x1;
            continue;
          }

          xu = x1;

          if ( s < m1 )
          {
            rv4[m1-1] = x1;
            continue;
          }
          rv4[s] = x1;
          if ( x1 < rv5[s-1] )
          {
            rv5[s-1] = x1;
          }
        }
/*
  K-th eigenvalue found.
*/
        rv5[k-1] = x1;
        k = k - 1;
        if ( k < m1 )
        {
          break;
        }
      }
    }
/*
  Order eigenvalues tagged with their submatrix associations.
*/
    s = r;
    r = r + m2 - m1 + 1;
    j = 1;
    k = m1;

    for ( l = 1; l <= r; l++ )
    {
      if ( j <= s )
      {
        if ( m2 < k )
        {
          break;
        }

        if ( w[l-1] <= rv5[k-1] )
        {
          j = j + 1;
          continue;
        }

        for ( i = l + s - j; l <= i; i-- )
        {
          w[i] = w[i-1];
          ind[i] = ind[i-1];
        }
      }

      w[l-1] = rv5[k-1];
      ind[l-1] = tag - 1;
      k = k + 1;
    }

    if ( n <= q )
    {
      break;
    }
  }

  free ( rv4 );
  free ( rv5 );

  return ierr;
}
/******************************************************************************/

int bqr ( int nm, int n, int mb, double a[], double *t, double *r )

/******************************************************************************/
/*
  Purpose:

    BQR is DUMMY CODE right now.
*/
{
  if ( true )
  {
    printf ( "\n" );
    printf ( "BQR - Fatal error!\n" );
    printf ( "  This is just DUMMY CODE right now.\n" );
    exit ( 1 );
  }
  return 0;
}
/******************************************************************************/

void cbabk2 ( int n, int low, int igh, double scale[], int m, double zr[], 
  double zi[] )

/******************************************************************************/
/*
  Purpose:

    CBABK2 finds eigenvectors by undoing the CBAL transformation.

  Discussion:

    CBABK2 forms the eigenvectors of a complex general
    matrix by back transforming those of the corresponding
    balanced matrix determined by CBAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, values determined by CBAL.

    Input, double SCALE[N], information determining the permutations
    and scaling factors used by CBAL.

    Input, int M, the number of eigenvectors to be back
    transformed.

    Input/output, double ZR[N*M], ZI[N*M].  On input, the real
    and imaginary parts, respectively, of the eigenvectors to be back
    transformed in their first M columns.  On output, the transformed
    eigenvectors.
*/
{
  int i;
  int ii;
  int j;
  int k;
  double s;

  if ( m == 0 )
  {
    return;
  }

  if ( igh != low )
  {
    for ( i = low; i <= igh; i++ )
    {
      s = scale[i];
      for ( j = 0; j < m; j++ )
      {
        zr[i+j*n] = zr[i+j*n] * s;
        zi[i+j*n] = zi[i+j*n] * s;
      }
    }
  }

  for ( ii = 0; ii < n; ii++ )
  {
    i = ii;
    if ( i < low || igh < i )
    {
      if ( i < low )
      {
        i = low - ii;
      }

      k = ( int ) scale[i];

      if ( k != i )
      {
        for ( j = 0; j < m; j++ )
        {
          s         = zr[i+j*n];
          zr[i+j*n] = zr[k+j*n];
          zr[k+j*n] = s;
          s         = zi[i+j*n];
          zi[i+j*n] = zi[k+j*n];
          zi[k+j*n] = s;
        }
      }
    }
  }
  return;
}
/******************************************************************************/

void cbal ( int n, double ar[], double ai[], int *low, int *igh, 
  double scale[] )

/******************************************************************************/
/*
  Purpose:

    CBAL balances a complex matrix before eigenvalue calculations.

  Discussion:

    CBAL balances a complex matrix and isolates eigenvalues whenever possible.

    Suppose that the principal submatrix in rows low through igh
    has been balanced, that P(J) denotes the index interchanged
    with J during the permutation step, and that the elements
    of the diagonal matrix used are denoted by D(I,J).  Then
      SCALE(J) = P(J),    for J = 1,...,LOW-1
               = D(J,J)       J = LOW,...,IGH
               = P(J)         J = IGH+1,...,N.
    The order in which the interchanges are made is N to IGH+1,
    then 1 to LOW-1.

    Note that 1 is returned for IGH if IGH is zero formally.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double AR[N*N], AI[N*N].  On input, the real and
    imaginary parts of the complex matrix to be balanced.  On output,
    the real and imaginary parts of the balanced matrix.

    Output, int *LOW, *IGH, are values such that AR(I,J)
    and AI(I,J) are zero if I is greater than J and either J=1,...,LOW-1 or
    I=IGH+1,...,N.

    Output, double SCALE[N], information determining the
    permutations and scaling factors used.
*/
{
  double b2;
  double c;
  double f;
  double g;
  int i;
  int j;
  bool jump;
  bool jump2;
  int k;
  int l;
  int m;
  bool noconv;
  double r;
  double radix;
  double s;
  double t;

  radix = 16.0;
  j = 0;
  m = 0;
  b2 = radix * radix;
  k = 0;
  l = n - 1;

  while ( true )
  {
    jump2 = true;

    for ( j = l; 0 <= j; j-- )
    {
      jump = true;
      for ( i = 0; i <= l; i++ )
      {
        if ( i != j )
        {
          if ( ar[j+i*n] != 0.0 || ai[j+i*n] != 0.0 )
          {
            jump = false;
            break;
          }
        }
      }

      if ( jump )
      {
        m = l;
        scale[m] = j;

        if ( j != m )
        {
          for ( i = 0; i <= l; i++ )
          {
            t         = ar[i+j*n];
            ar[i+j*n] = ar[i+m*n];
            ar[i+m*n] = t;
            t         = ai[i+j*n];
            ai[i+j*n] = ai[i+m*n];
            ai[i+m*n] = t;
          }
          for ( i = k; i < n; i++ )
          {
            t         = ar[j+i*n];
            ar[j+i*n] = ar[m+i*n];
            ar[m+i*n] = t;
            t         = ai[j+i*n];
            ai[j+i*n] = ai[m+i*n];
            ai[m+i*n] = t;
          }
        }

        if ( l == 0 )
        {
          *low = k;
          *igh = l;
          return;
        }

        l = l - 1;
        jump2 = false;
        break;
      }
    }

    if ( jump2 )
    {
      break;
    }
  }
/*
  Search for columns isolating an eigenvalue and push them left.
*/
  while ( true )
  {
    jump2 = true;

    for ( j = k; j <= l; j++ )
    {
      jump = true;
      for ( i = k; i <= l; i++ )
      {
        if ( i != j )
        {
          if ( ar[i+j*n] != 0.0 || ai[i+j*n] != 0.0 )
          {
            jump = false;
            break;
          }
        }
      }

      if ( jump )
      {
        m = k;
        scale[m] = j;

        if ( j != m )
        {
          for ( i = 0; i <= l; i++ )
          {
            t         = ar[i+j*n];
            ar[i+j*n] = ar[i+m*n];
            ar[i+m*n] = t;
            t         = ai[i+j*n];
            ai[i+j*n] = ai[i+m*n];
            ai[i+m*n] = t;
          }
          for ( i = k; i < n; i++ )
          {
            t         = ar[j+i*n];
            ar[j+i*n] = ar[m+i*n];
            ar[m+i*n] = t;
            t         = ai[j+i*n];
            ai[j+i*n] = ai[m+i*n];
            ai[m+i*n] = t;
          }
        }
        k = k + 1;
        jump2 = false;
        break;
      }
    }

    if ( jump2 )
    {
      break;
    }
  }
/*
  Now balance the submatrix in rows k to l.
*/
  for ( i = k; i <= l; i++ )
  {
    scale[i] = 1.0;
  }
/*
  Iterative loop for norm reduction.
*/
  while ( true )
  {
    noconv = false;
    for ( i = k; i <= l; i++ )
    {
      c = 0.0;
      r = 0.0;
      for ( j = k; j <= l; j++ )
      {
        if ( j != i )
        {
          c = c + fabs ( ar[j+i*n] ) + fabs ( ai[j+i*n] );
          r = r + fabs ( ar[i+j*n] ) + fabs ( ai[i+j*n] );
        }
      }
/*
  Guard against zero C or R due to underflow.
*/
      if ( c != 0.0 && r != 0.0 )
      {
        g = r / radix;
        f = 1.0;
        s = c + r;
        while ( c < g )
        {
          f = f * radix;
          c = c * b2;
        }
        g = r * radix;
        while  ( g <= c )
        {
          f = f / radix;
          c = c / b2;
        }
/*
  Now balance.
*/
        if ( ( c + r ) / f < 0.95 * s )
        {
          g = 1.0 / f;
          scale[i] = scale[i] * f;
          noconv = true;
          for ( j = k; j < n; j++ )
          {
            ar[i+j*n] = ar[i+j*n] * g;
            ai[i+j*n] = ai[i+j*n] * g;
          }
          for ( j = 0; j <= l; j++ )
          {
            ar[j+i*n] = ar[j+i*n] * f;
            ai[j+i*n] = ai[j+i*n] * f;
          }
        }
      }
    }

    if ( ! noconv )
    {
      break;
    }
  }

  *low = k;
  *igh = l;

  return;
}
/******************************************************************************/

void cdiv ( double ar, double ai, double br, double bi, double *cr, double *ci )

/******************************************************************************/
/*
  Purpose:

    CDIV emulates complex division, using real arithmetic.

  Discussion:

    CDIV performs complex division:

      (CR,CI) = (AR,AI) / (BR,BI)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2009

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, double AR, AI, the real and imaginary parts of
    the numerator.

    Input, double BR, BI, the real and imaginary parts of
    the denominator.

    Output, double *CR, *CI, the real and imaginary parts of
    the result.
*/
{
  double ais;
  double ars;
  double bis;
  double brs;
  double s;

  s = fabs ( br ) + fabs ( bi );

  ars = ar / s;
  ais = ai / s;
  brs = br / s;
  bis = bi / s;

  s = brs * brs + bis * bis;
  *cr = ( ars * brs + ais * bis ) / s;
  *ci = ( ais * brs - ars * bis ) / s;

  return;
}
/******************************************************************************/

int cg_lr ( int n, double ar[], double ai[], double wr[], double wi[], 
  bool matz, double zr[], double zi[] )

/******************************************************************************/
/*
  Purpose:

    CG_LR gets eigenvalues and eigenvectors of a complex general matrix.

  Discussion:

    CG_LR calls EISPACK routines to find the eigenvalues and eigenvectors 
    of a complex general matrix, using elementary transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 March 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double AR[N,N], AI[N,N].  On input, the real and
    imaginary parts of the complex matrix.  On output, AR and AI
    have been overwritten by other information.

    Output, double WR[N], WI[N], the real and imaginary parts
    of the eigenvalues.

    Input, bool MATZ, is false if only eigenvalues are desired, and
    true if both eigenvalues and eigenvectors are to be computed.

    Output, double ZR[N*N], ZI[N*N], the real and imaginary parts,
    respectively, of the eigenvectors, if MATZ is true.

    Output, int CG_LR, an error completion code described in 
    the documentation for COMLR and COMLR2.  The normal completion code 
    is zero.
*/
{
  double *fv1;
  int i;
  int ierr;
  int *inter;
  int is1;
  int is2;

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );

  cbal ( n, ar, ai, &is1, &is2, fv1 );

  inter = ( int * ) malloc ( n * sizeof ( int ) );
  for ( i = 0; i < n; i++ )
  {
    inter[i] = -1;
  }

  comhes ( n, is1, is2, ar, ai, inter );

  if ( ! matz )
  {
    ierr = comlr ( n, is1, is2, ar, ai, wr, wi );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "CG_LR - Fatal error!\n" );
      printf ( "  Nonzero error return from COMQR.\n" );
      free ( fv1 );
      free ( inter );
      return ierr;
    }
  }
  else
  {
    ierr = comlr2 ( n, is1, is2, inter, ar, ai, wr, wi, zr, zi );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "CG_LR - Fatal error!\n" );
      printf ( "  Nonzero error return from COMQR2.\n" );
      free ( fv1 );
      free ( inter );
      return ierr;
    }

    cbabk2 ( n, is1, is2, fv1, n, zr, zi );
  }

  free ( fv1 );
  free ( inter );

  return ierr;
}
/******************************************************************************/

int cg_qr ( int n, double ar[], double ai[], double wr[], double wi[], 
  bool matz, double zr[], double zi[] )

/******************************************************************************/
/*
  Purpose:

    CG_QR gets eigenvalues and eigenvectors of a complex general matrix.

  Discussion:

    CG_QR calls EISPACK routines to find the eigenvalues and eigenvectors 
    of a complex general matrix, using unitary transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double AR[N,N], AI[N,N].  On input, the real and
    imaginary parts of the complex matrix.  On output, AR and AI
    have been overwritten by other information.

    Output, double WR[N], WI[N], the real and imaginary parts
    of the eigenvalues.

    Input, bool MATZ, is false if only eigenvalues are desired, and
    true if both eigenvalues and eigenvectors are to be computed.

    Output, double ZR[N,N], ZI[N,N], the real and imaginary parts,
    respectively, of the eigenvectors, if MATZ is true.

    Output, int CG_QR, an error completion code described in 
    the documentation for COMQR and COMQR2.  The normal completion code 
    is zero.
*/
{
  double *fv1;
  double *fv2;
  double *fv3;
  int ierr;
  int is1;
  int is2;

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );

  cbal ( n, ar, ai, &is1, &is2, fv1 );

  fv2 = ( double * ) malloc ( n * sizeof ( double ) );
  fv3 = ( double * ) malloc ( n * sizeof ( double ) );

  corth ( n, is1, is2, ar, ai, fv2, fv3 );

  if ( ! matz )
  {
    ierr = comqr ( n, is1, is2, ar, ai, wr, wi );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "CG_QR - Fatal error!\n" );
      printf ( "  Nonzero error return from COMQR.\n" );
      free ( fv1 );
      free ( fv2 );
      free ( fv3 );
      return ierr;
    }
  }
  else
  {
    ierr = comqr2 ( n, is1, is2, fv2, fv3, ar, ai, wr, wi, zr, zi );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "CG_QR - Fatal error!\n" );
      printf ( "  Nonzero error return from COMQR2.\n" );
      free ( fv1 );
      free ( fv2 );
      free ( fv3 );
      return ierr;
    }

    cbabk2 ( n, is1, is2, fv1, n, zr, zi );
  }

  free ( fv1 );
  free ( fv2 );
  free ( fv3 );

  return ierr;
}
/******************************************************************************/

int ch ( int n, double ar[], double ai[], double w[], bool matz, double zr[], 
  double zi[] )

/******************************************************************************/
/*
  Purpose:

    CH gets eigenvalues and eigenvectors of a complex Hermitian matrix.

  Discussion:

    CH calls the recommended sequence of routines from the
    EISPACK eigensystem package to find the eigenvalues and eigenvectors
    of a complex hermitian matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double AR[N,N], AI[N,N].  On input, the real and
    imaginary parts of the complex matrix.  On output, AR and AI
    have been overwritten by other information.

    Output, double W[N], the eigenvalues in ascending order.

    Input, bool MATZ, is false if only eigenvalues are desired, and
    true if both eigenvalues and eigenvectors are to be computed.

    Output, double ZR[N,N], ZI[N,N], the real and imaginary parts,
    respectively, of the eigenvectors, if MATZ is true.

    Output, int CH, an error completion code described in 
    the documentation for TQLRAT and TQL2.  The normal completion code is zero.
*/
{
  double *fm1;
  double *fv1;
  double *fv2;
  int ierr;

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );
  fv2 = ( double * ) malloc ( n * sizeof ( double ) );
  fm1 = ( double * ) malloc ( 2 * n * sizeof ( double ) );

  htridi ( n, ar, ai, w, fv1, fv2, fm1 );

  if ( ! matz )
  {
    ierr = tqlrat ( n, w, fv2 );
    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );
      free ( fm1 );
      return ierr;
    }
  }
  else
  {
    r8mat_identity ( n, zr );

    ierr = tql2 ( n, w, fv1, zr );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );
      free ( fm1 );
      return ierr;
    }

    htribk ( n, ar, ai, fm1, n, zr, zi );
  }

  free ( fv1 );
  free ( fv2 );
  free ( fm1 );

  return ierr;
}
/******************************************************************************/

int ch3 ( int n, double a[], double w[], bool matz, double zr[], double zi[] )

/******************************************************************************/
/*
  Purpose:

    CH3 gets eigenvalues and eigenvectors of a complex Hermitian matrix.

  Discussion:

    CH3 calls the recommended sequence of routines from the
    EISPACK eigensystem package to find the eigenvalues and eigenvectors
    of a complex hermitian matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 February 2018

  Author:

    John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double A[N,N].  On input, the lower triangle of 
    the complex hermitian input matrix.  The real parts of the matrix elements 
    are stored in the full lower triangle of A, and the imaginary parts are 
    stored in the transposed positions of the strict upper triangle of A.  No 
    storage is required for the zero imaginary parts of the diagonal elements.
    On output, A contains information about the unitary transformations
    used in the reduction.

    Output, double W[N], the eigenvalues in ascending order.

    Input, bool MATZ, is false if only eigenvalues are desired, and
    true if both eigenvalues and eigenvectors are to be computed.

    Output, double ZR[N,N], ZI[N,N], the real and imaginary parts,
    respectively, of the eigenvectors, if MATZ is true.

    Output, int CH, an error completion code described in 
    the documentation for TQLRAT and TQL2.  The normal completion code is zero.
*/
{
  double *fm1;
  double *fv1;
  double *fv2;
  int ierr;

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );
  fv2 = ( double * ) malloc ( n * sizeof ( double ) );
  fm1 = ( double * ) malloc ( 2 * n * sizeof ( double ) );

  htrid3 ( n, a, w, fv1, fv2, fm1 );

  if ( ! matz )
  {
    ierr = tqlrat ( n, w, fv2 );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );
      free ( fm1 );
      return ierr;
    }
  }
  else
  {
    r8mat_identity ( n, zr );
    r8mat_zeros ( n, n, zi );

    ierr = tql2 ( n, w, fv1, zr );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );
      free ( fm1 );
      return ierr;
    }

    htrib3 ( n, a, fm1, n, zr, zi );
  }

  free ( fv1 );
  free ( fv2 );
  free ( fm1 );

  return ierr;
}
/******************************************************************************/

int cinvit ( int n, double ar[], double ai[], double wr[], double wi[], 
  bool select[], int mm, int *m, double zr[], double zi[] )

/******************************************************************************/
/*
  Purpose:

    CINVIT gets eigenvectors from eigenvalues for a complex Hessenberg matrix.

  Discussion:

    CINVIT finds those eigenvectors of a complex upper Hessenberg matrix 
    corresponding to specified eigenvalues, using inverse iteration.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 March 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double AR[N,N], AI[N,N], the real and imaginary parts of
    the complex Hessenberg matrix.

    Input/output, double WR[N], WI[N].  On input, the real and
    imaginary parts of the eigenvalues of the matrix.  The eigenvalues must
    be stored in a manner identical to that of COMLR, which
    recognizes possible splitting of the matrix.  On output, WR may have been
    altered since close eigenvalues are perturbed slightly in searching for
    independent eigenvectors.

    Input, logical SELECT[N], specifies the eigenvectors to be found.  The
    eigenvector corresponding to the J-th eigenvalue is specified by
    setting SELECT(J) to TRUE.

    Input, int MM, an upper bound for the number of 
    eigenvectors to be found.

    Output, int M, the number of eigenvectors actually found.

    Output, double ZR[N,MM], ZI[N,MM], the real and imaginary parts
    of the eigenvectors.  The eigenvectors are normalized so that the
    component of largest magnitude is 1.
    Any vector which fails the acceptance test is set to zero.

    Output, int CINVIT, error flag.
    0, for normal return,
    -(2*N+1), if more than MM eigenvectors have been specified,
    -K, if the iteration corresponding to the K-th value fails,
    -(N+K), if both error situations occur.

*/
{
  double eps3 = 0.0 ;
  double growto = 0.0 ;
  int i;
  int ierr;
  double ilambd;
  int its;
  int j = 0 ;
  int k;
  int mp;
  double norm;
  double normv;
  bool repeat;
  double rlambd;
  double *rm1;
  double *rm2;
  double *rv1;
  double *rv2;
  int s;
  double t;
  double ti;
  double tr;
  int uk;
  double ukroot = 0.0 ;
  double x;
  double y;

  rm1 = ( double * ) malloc ( n * n * sizeof ( double ) );
  rm2 = ( double * ) malloc ( n * n * sizeof ( double ) );
  rv1 = ( double * ) malloc ( n * sizeof ( double ) );
  rv2 = ( double * ) malloc ( n * sizeof ( double ) );

  ierr = 0;
  uk = 0;
  s = 1;

  for ( k = 0; k < n; k++ )
  {
    if ( ! select[k] )
    {
      continue;
    }

    if ( mm < s )
    {
      if ( ierr != 0 )
      {
        ierr = ierr - n;
      }

      if ( ierr == 0 )
      {
        ierr = - ( 2 * n + 1 );
      }
      *m = s - 1;
      free ( rm1 );
      free ( rm2 );
      free ( rv1 );
      free ( rv2 );
      return ierr;
    }

    if ( uk < k )
    {
/*
  Check for possible splitting.
*/
      for ( uk = k; uk < n; uk++ )
      {
        if ( uk == n - 1 )
        {
          break;
        }
        if ( ar[uk+1+uk*n] == 0.0 && ai[uk+1+uk*n] == 0.0 )
        {
          break;
        }
      }
/*
  Compute infinity norm of leading UK by UK (Hessenberg) matrix.
*/
      norm = 0.0;
      for ( i = 0; i <= uk; i++ )
      {
        x = 0.0;
        for ( j = i4_max ( 0, i - 1 ); j <= uk; j++ )
        {
          x = x + pythag ( ar[i+j*n], ai[i+j*n] );
        }
        norm = r8_max ( norm, x );
      }
/*
  EPS3 replaces zero pivot in decomposition
  and close roots are modified by EPS3.
*/
      if ( norm == 0.0 )
      {
        norm = 1.0;
      }

      eps3 = fabs ( norm ) * r8_epsilon ( );
/*
  GROWTO is the criterion for growth.
*/
      ukroot = ( double ) ( uk );
      ukroot = sqrt ( ukroot );
      growto = 0.1 / ukroot;
    }

    rlambd = wr[k];
    ilambd = wi[k];
/*
  Perturb eigenvalue if it is close to any previous eigenvalue.
*/
    if ( 0 < k )
    {
      while ( true )
      {
        repeat = false;

        for ( i = k - 1; 0 <= i; i-- )
        {
          if ( select[i] &&
               fabs ( wr[i] - rlambd ) < eps3 &&
               fabs ( wi[i] - ilambd ) < eps3 )
          {
            rlambd = rlambd + eps3;
            repeat = true;
            break;
          }
        }

        if ( ! repeat )
        {
          wr[k] = rlambd;
          break;
        }
      }
    }

    mp = 0;

    for ( i = 0; i <= uk; i++ )
    {
      for ( j = i4_max ( 0, i - 1 ); j <= uk; j++ )
      {
        rm1[i+j*n] = ar[i+j*n];
        rm2[i+j*n] = ai[i+j*n];
      }
      rm1[i+i*n] = rm1[i+i*n] - rlambd;
      rm2[i+i*n] = rm2[i+i*n] - ilambd;
      rv1[i] = eps3;
    }
/*
  Triangular decomposition with interchanges, replacing zero pivots by eps3.
*/
    for ( i = 1; i <= uk; i++ )
    {
      mp = i - 1;

      if ( pythag ( rm1[mp+mp*n], rm2[mp+mp*n] ) < 
           pythag ( rm1[i+mp*n], rm2[i+mp*n] ) )
      {
        for ( j = i - 1; j <= uk; j++ )
        {
          t           = rm1[i+j*n];
          rm1[i+j*n]  = rm1[mp+j*n];
          rm1[mp+j*n] = t;
          t           = rm2[i+j*n];
          rm2[i+j*n]  = rm2[mp+j*n];
          rm2[mp+j*n] = t;
        }
      }

      if ( rm1[mp+mp*n] == 0.0 && rm2[mp+mp*n] == 0.0 )
      {
        rm1[mp+mp*n] = eps3;
      }

      cdiv ( rm1[i+mp*n], rm2[i+mp*n], rm1[mp+mp*n], rm2[mp+mp*n], &x, &y );

      if ( x != 0.0 || y != 0.0 )
      {
        for ( j = i; j <= uk; j++ )
        {
          rm1[i+j*n] = rm1[i+j*n] - x * rm1[mp+j*n] + y * rm2[mp+j*n];
          rm2[i+j*n] = rm2[i+j*n] - x * rm2[mp+j*n] - y * rm1[mp+j*n];
        }
      }
    }

    if ( rm1[uk+uk*n] == 0.0 && rm2[uk+uk*n] == 0.0 )
    {
      rm1[uk+uk*n] = eps3;
    }

    its = 0;
/*
  Back substitution.
*/
    while ( true )
    {
      for ( i = uk; 0 <= i; i-- )
      {
        x = rv1[i];
        y = 0.0;
        for ( j = i + 1; j <= uk; j++ )
        {
          x = x - rm1[i+j*n] * rv1[j] + rm2[i+j*n] * rv2[j];
          y = y - rm1[i+j*n] * rv2[j] - rm2[i+j*n] * rv1[j];
        }
        cdiv ( x, y, rm1[i+i*n], rm2[i+i*n], &tr, &ti );
        rv1[i] = tr;
        rv2[i] = ti;
      }
/*
  Acceptance test for eigenvector and normalization.
*/
      its = its + 1;
      norm = 0.0;
      normv = 0.0;

      for ( i = 0; i <= uk; i++ )
      {
        x = pythag ( rv1[i], rv2[i] );
        if ( normv < x )
        {
          normv = x;
          j = i;
        }
        norm = norm + x;
      }

      if ( growto <= norm )
      {
/*
  Accept vector.
*/
        x = rv1[j];
        y = rv2[j];

        for ( i = 0; i <= uk; i++ )
        {
          cdiv ( rv1[i], rv2[i], x, y, &tr, &ti );
          zr[i+(s-1)*n] = tr;
          zi[i+(s-1)*n] = ti;
        }

        if ( uk != n - 1 )
        {
          j = uk + 1;
          for ( i = j; i < n; i++ )
          {
            zr[i+(s-1)*n] = 0.0;
            zi[i+(s-1)*n] = 0.0;
          }
        }
        s = s + 1;
        break;
      }
/*
  Choose a new starting vector.
*/
      else if ( its < uk )
      {
        x = ukroot;
        y = eps3 / ( x + 1.0 );

        rv1[0] = eps3;
        for ( i = 1; i <= uk; i++ )
        {
          rv1[i] = y;
        }
        j = uk - its;
        rv1[j] = rv1[j] - eps3 * x;
      }
/*
  Error: unaccepted eigenvector.
*/
      else
      {
        j = 0;
        ierr = - k;
/*
  Set remaining vector components to zero.
*/
        for ( i = j; i < n; i++ )
        {
          zr[i+(s-1)*n] = 0.0;
          zi[i+(s-1)*n] = 0.0;
        }
        s = s + 1;
        break;
      }
    }
  }

  *m = s - 1;

  free ( rm1 );
  free ( rm2 );
  free ( rv1 );
  free ( rv2 );

  return ierr;
}
/******************************************************************************/

void combak ( int n, int low, int igh, double ar[], double ai[], int inter[], 
  int m, double zr[], double zi[] )

/******************************************************************************/
/*
  Purpose:

    COMBAK determines eigenvectors by undoing the COMHES transformation.

  Discussion:

    COMBAK forms the eigenvectors of a complex general
    matrix by back transforming those of the corresponding
    upper Hessenberg matrix determined by COMHES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing
    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = to the order
    of the matrix.

    Input, double AR[N,IGH], AI[N,IGH], the multipliers which
    were used in the reduction by COMHES in their lower triangles below
    the subdiagonal.

    Input, int INTER[IGH], information on the rows and
    columns interchanged in the reduction by COMHES.

    Input, int M, the number of eigenvectors to be back
    transformed.

    Input/output, double ZR[N,M], ZI[N,M].  On input, the real
    and imaginary parts of the eigenvectors to be back transformed.  On
    output, the real and imaginary parts of the transformed eigenvectors.
*/
{
  int i;
  int j;
  int mm;
  int mp;
  double t;
  double xi;
  double xr;

  if ( m == 0 )
  {
    return;
  }

  if ( igh - 1 < low + 1 )
  {
    return;
  }

  for ( mm = low + 1; mm <= igh - 1; mm++ )
  {
     mp = low + igh - mm;

     for ( i = mp + 1; i <= igh; i++ )
     {
        xr = ar[i+(mp-1)*n];
        xi = ai[i+(mp-1)*n];

        if ( xr != 0.0 || xi != 0.0 )
        {
          for ( j = 0; j < m; j++ )
          {
            zr[i+j*n] = zr[i+j*n] + xr * zr[mp+j*n] - xi * zi[mp+j*n];
            zi[i+j*n] = zi[i+j*n] + xr * zi[mp+j*n] + xi * zr[mp+j*n];
          }
       }
     }

     i = inter[mp];

     if ( i != mp )
     {
       for ( j = 0; j < m; j++ )
       {
         t          = zr[i+j*n];
         zr[i+j*n]  = zr[mp+j*n];
         zr[mp+j*n] = t;

         t          = zi[i+j*n];
         zi[i+j*n]  = zi[mp+j*n];
         zi[mp+j*n] = t;
       }
     }
  }

  return;
}
/******************************************************************************/

void comhes ( int n, int low, int igh, double ar[], double ai[], int inter[] )

/******************************************************************************/
/*
  Purpose:

    COMHES transforms a complex general matrix to upper Hessenberg form.

  Discussion:

    COMHES is given a complex general matrix and reduces a submatrix in rows 
    and columns LOW through IGH to upper Hessenberg form by
    stabilized elementary similarity transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing
    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.

    Input/output, double AR[N,N], AI[N,N].  On input, the real and
    imaginary parts of the complex input matrix.  On output, the real and
    imaginary parts of the Hessenberg matrix.  The multipliers which were
    used in the reduction are stored in the remaining triangles under the
    Hessenberg matrix.

    Output, int INTER[IGH], information on the rows and 
    columns interchanged in the reduction.
*/
{
  int i;
  int j;
  int m;
  double t;
  double xi;
  double xr;
  double yi;
  double yr;

  for ( m = low + 1; m <= igh - 1; m++ )
  {
/*
  Choose the pivot I.
*/
    xr = 0.0;
    xi = 0.0;
    i = m;

    for ( j = m; j <= igh; j++ )
    {
      if ( fabs ( xr )            + fabs ( xi ) < 
           fabs ( ar[j+(m-1)*n] ) + fabs ( ai[j+(m-1)*n] ) )
      {
        xr = ar[j+(m-1)*n];
        xi = ai[j+(m-1)*n];
        i = j;
      }
    }

    inter[m] = i;
/*
  Interchange rows and columns of AR and AI.
*/
    if ( i != m )
    {
      for ( j = m - 1; j < n; j++ )
      {
        t         = ar[i+j*n];
        ar[i+j*n] = ar[m+j*n];
        ar[m+j*n] = t;
        t         = ai[i+j*n];
        ai[i+j*n] = ai[m+j*n];
        ai[m+j*n] = t;
      }

      for ( j = 0; j <= igh; j++ )
      {
        t         = ar[j+i*n];
        ar[j+i*n] = ar[j+m*n];
        ar[j+m*n] = t;
        t         = ai[j+i*n];
        ai[j+i*n] = ai[j+m*n];
        ai[j+m*n] = t;
      }
    }
/*
  Carry out the transformation.
*/
    if ( xr != 0.0 || xi != 0.0 )
    {
      for ( i = m + 1; i <= igh; i++ )
      {
        yr = ar[i+(m-1)*n];
        yi = ai[i+(m-1)*n];

        if ( yr != 0.0 || yi != 0.0 )
        {
          cdiv ( yr, yi, xr, xi, &yr, &yi );
          ar[i+(m-1)*n] = yr;
          ai[i+(m-1)*n] = yi;

          for ( j = m; j < n; j++ )
          {
            ar[i+j*n] = ar[i+j*n] - yr * ar[m+j*n] + yi * ai[m+j*n];
            ai[i+j*n] = ai[i+j*n] - yr * ai[m+j*n] - yi * ar[m+j*n];
          }

          for ( j = 0; j <= igh; j++ )
          {
            ar[j+m*n] = ar[j+m*n] + yr * ar[j+i*n] - yi * ai[j+i*n];
            ai[j+m*n] = ai[j+m*n] + yr * ai[j+i*n] + yi * ar[j+i*n];
          }
        }
      }
    }
  }

  return;
}
/******************************************************************************/

int comlr ( int n, int low, int igh, double hr[], double hi[], double wr[], 
  double wi[] )

/******************************************************************************/
/*
  Purpose:

    COMLR gets all eigenvalues of a complex upper Hessenberg matrix.

  Discussion:

    COMLR finds the eigenvalues of a complex upper Hessenberg matrix by the 
    modified LR method.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing
    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.

    Input/output, double HR[N,N], HI[N,N].  On input, the real and
    imaginary parts of the complex upper Hessenberg matrix.  Their lower
    triangles below the subdiagonal contain the multipliers which were used
    in the reduction by COMHES if performed.  On output, the upper Hessenberg
    portions of HR and HI have been destroyed.  Therefore, they must be
    saved before calling COMLR if subsequent calculation of eigenvectors
    is to be performed.

    Output, double WR[N], WI[N], the real and imaginary parts of the
    eigenvalues.  If an error break; is made, the eigenvalues should be correct
    for indices IERR+1,...,N.

    Output, int COMLR, error flag.
    0, for normal return,
    J, if the limit of 30*N iterations is exhausted while the J-th
      eigenvalue is being sought.
*/
{
  double ai;
  double ar;
  int en;
  int i;
  int ierr;
  int itn;
  int its;
  int j;
  int l;
  int m;
  double si;
  double sr;
  double t;
  double ti;
  double tr;
  double tst1;
  double tst2;
  double xi;
  double xr;
  double yi;
  double yr;
  double zzi;
  double zzr;

  ierr = 0;
/*
  Store roots isolated by CBAL.
*/
  for ( i = 0; i < n; i++ )
  {
    if ( i < low || igh < i )
    {
      wr[i] = hr[i+i*n];
      wi[i] = hi[i+i*n];
    }
  }

  en = igh;
  tr = 0.0;
  ti = 0.0;
  itn = 30 * n;
/*
  Search for next eigenvalue.
*/
  if ( en < low )
  {
    return ierr;
  }

  its = 0;
/*
  Look for single small sub-diagonal element.
*/
  while ( true )
  {
    for ( l = en; low <= l; l-- )
    {
      if ( l == low )
      {
        break;
      }

      tst1 = fabs ( hr[l-1+(l-1)*n] ) + fabs ( hi[l-1+(l-1)*n] ) + fabs ( hr[l+l*n] ) 
        + fabs ( hi[l+l*n] );
      tst2 = tst1 + fabs ( hr[l+(l-1)*n] ) + fabs ( hi[l+(l-1)*n] );

      if ( tst2 == tst1 )
      {
        break;
      }
    }
/*
  A root found.
*/
    if ( l == en )
    {
      wr[en] = hr[en+en*n] + tr;
      wi[en] = hi[en+en*n] + ti;
      en = en - 1;
      if ( en < low )
      {
        return ierr;
      }
      its = 0;
      continue;
    }

    if ( itn == 0 )
    {
      ierr = en;
      return ierr;
    }

    if ( its == 10 || its == 20 )
    {
      sr = fabs ( hr[en+(en-1)*n] ) + fabs ( hr[en-1+(en-2)*n] );
      si = fabs ( hi[en+(en-1)*n] ) + fabs ( hi[en-1+(en-2)*n] );
    }
    else
    {
      sr = hr[en+en*n];
      si = hi[en+en*n];
      xr = hr[en-1+en*n] * hr[en+(en-1)*n] - hi[en-1+en*n] * hi[en+(en-1)*n];
      xi = hr[en-1+en*n] * hi[en+(en-1)*n] + hi[en-1+en*n] * hr[en+(en-1)*n];

      if ( xr != 0.0 || xi != 0.0 )
      {
        yr = ( hr[en-1+(en-1)*n] - sr) / 2.0;
        yi = ( hi[en-1+(en-1)*n] - si) / 2.0;
        ar = yr * yr - yi * yi + xr;
        ai = 2.0 * yr * yi + xi;
        csroot ( ar, ai, &zzr, &zzi );

        if ( yr * zzr + yi * zzi < 0.0 )
        {
          zzr = - zzr;
          zzi = - zzi;
        }

        ar = yr + zzr;
        ai = yi + zzi;
        cdiv ( xr, xi, ar, ai, &xr, &xi );
        sr = sr - xr;
        si = si - xi;
      }
    }

    for ( i = low; i <= en; i++ )
    {
      hr[i+i*n] = hr[i+i*n] - sr;
      hi[i+i*n] = hi[i+i*n] - si;
    }

    tr = tr + sr;
    ti = ti + si;
    its = its + 1;
    itn = itn - 1;
/*
  Look for two consecutive small sub-diagonal elements.
*/
    xr = fabs ( hr[en-1+(en-1)*n] ) + fabs ( hi[en-1+(en-1)*n] );
    yr = fabs ( hr[en+(en-1)*n] ) + fabs ( hi[en+(en-1)*n] );
    zzr = fabs ( hr[en+en*n] ) + fabs ( hi[en+en*n] );

    for ( m = en - 1; l <= m; m-- )
    {
      if ( m == l )
      {
        break;
      }

      yi = yr;
      yr = fabs ( hr[m+(m-1)*n] ) + fabs ( hi[m+(m-1)*n] );
      xi = zzr;
      zzr = xr;
      xr = fabs ( hr[m-1+(m-1)*n] ) + fabs ( hi[m-1+(m-1)*n] );
      tst1 = zzr / yi * ( zzr + xr + xi );
      tst2 = tst1 + yr;
      if ( tst2 == tst1 ) 
      {
        break;
      }
    }
/*
  Triangular decomposition H=L*R.
*/
    for ( i = m + 1; i <= en; i++ )
    {
      xr = hr[i-1+(i-1)*n];
      xi = hi[i-1+(i-1)*n];
      yr = hr[i+(i-1)*n];
      yi = hi[i+(i-1)*n];

      if (  fabs ( yr ) + fabs ( yi ) <= fabs ( xr ) + fabs ( xi ) )
      {
        cdiv ( yr, yi, xr, xi, &zzr, &zzi );
        wr[i] = - 1.0;
      }
/*
  Interchange rows of HR and HI.
*/
      else
      {
        for ( j = i - 1; j <= en; j++ )
        {
          t           = hr[i-1+j*n];
          hr[i-1+j*n] = hr[i+j*n];
          hr[i+j*n]   = t;
          t           = hi[i-1+j*n];
          hi[i-1+j*n] = hi[i+j*n];
          hi[i+j*n]   = t;
        }

        cdiv ( xr, xi, yr, yi, &zzr, &zzi );
        wr[i] = 1.0;
      }

      hr[i+(i-1)*n] = zzr;
      hi[i+(i-1)*n] = zzi;

      for ( j = i; j <= en; j++ )
      {
        hr[i+j*n] = hr[i+j*n] - zzr * hr[i-1+j*n] + zzi * hi[i-1+j*n];
        hi[i+j*n] = hi[i+j*n] - zzr * hi[i-1+j*n] - zzi * hr[i-1+j*n];
      }
    }
/*
  Composition R*L=H.
*/
    for ( j = m + 1; j <= en; j++ )
    {
      xr = hr[j+(j-1)*n];
      xi = hi[j+(j-1)*n];
      hr[j+(j-1)*n] = 0.0;
      hi[j+(j-1)*n] = 0.0;
/*
  Interchange columns of HR and HI, if necessary.
*/
      if ( 0.0 < wr[j] )
      {
        for ( i = l; i <= j; i++ )
        {
          t             = hr[i+(j-1)*n];
          hr[i+(j-1)*n] = hr[i+j*n];
          hr[i+j*n]     = t;
          t             = hi[i+(j-1)*n];
          hi[i+(j-1)*n] = hi[i+j*n];
          hi[i+j*n]     = t;
        }
      }

      for ( i = l; i <= j; i++ )
      {
        hr[i+(j-1)*n] = hr[i+(j-1)*n] + xr * hr[i+j*n] - xi * hi[i+j*n];
        hi[i+(j-1)*n] = hi[i+(j-1)*n] + xr * hi[i+j*n] + xi * hr[i+j*n];
      }
    }
  }

  return ierr;
}
/******************************************************************************/

int comlr2 ( int n, int low, int igh, int inter[], double hr[], double hi[], 
  double wr[], double wi[], double zr[], double zi[] )

/******************************************************************************/
/*
  Purpose:

    COMLR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.

  Discussion:

    COMLR2 finds the eigenvalues and eigenvectors of a complex
    upper Hessenberg matrix by the modified LR method.  

    The eigenvectors of a complex general matrix can also be found if 
    COMHES has been used to reduce this general matrix to Hessenberg form.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing
    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.

    Input, int INTER[IGH], information on the rows and columns
    interchanged in the reduction by COMHES, if performed.  If the
    eigenvectors of the Hessenberg matrix are desired, set INTER(J)=J for these
    elements.

    Input/output, double HR[N,N], HI[N,N].  On input, the real
    and imaginary parts of the complex upper Hessenberg matrix.  Their lower
    triangles below the subdiagonal contain the multipliers which were used in
    the reduction by COMHES, if performed.  If the eigenvectors of the
    Hessenberg matrix are desired, these elements must be set to zero.  On
    output, the upper Hessenberg portions of HR and HI have been destroyed,
    but the location HR(1,1) contains the norm of the triangularized matrix.

    Output, double WR[N], WI[N], the real and imaginary parts of the
    eigenvalues.  If an error exit is made, the eigenvalues should be
    correct for indices IERR+1,...,N.

    Output, double ZR[N,N], ZI[N,N], the real and imaginary parts
    of the eigenvectors.  The eigenvectors are unnormalized.  If an error exit
    is made, none of the eigenvectors has been found.

    Output, int COMLR2, error flag.
    0, for normal return,
    J, if the limit of 30*N iterations is exhausted while the J-th
      eigenvalue is being sought.
*/
{
  double ai;
  double ar;
  int en;
  int enm1;
  int i;
  int ierr;
  int itn;
  int its;
  int j;
  int k;
  int l;
  int m;
  double norm;
  double si;
  double sr;
  double t;
  double ti;
  double tr;
  double tst1;
  double tst2;
  double xi;
  double xr;
  double yi;
  double yr;
  double zzi;
  double zzr;

  ierr = 0;
/*
  Initialize the eigenvector matrix.
*/
  r8mat_identity ( n, zr );
  r8mat_zeros ( n, n, zi );
/*
  Form the matrix of accumulated transformations from the information left
  by COMHES.
*/
  for ( i = igh - 1; low + 1 <= i; i-- )
  {
    for ( k = i + 1; k <= igh; k++ )
    {
      zr[k+i*n] = hr[k+(i-1)*n];
      zi[k+i*n] = hi[k+(i-1)*n];
    }

    j = inter[i];

    if ( i != j )
    {
      for ( k = i; k <= igh; k++ )
      {
        zr[i+k*n] = zr[j+k*n];
        zi[i+k*n] = zi[j+k*n];
        zr[j+k*n] = 0.0;
        zi[j+k*n] = 0.0;
      }
      zr[j+i*n] = 1.0;
    }
  }
/*
  Store roots isolated by CBAL.
*/
  for ( i = 0; i < n; i++ )
  {
    if ( i < low || igh < i )
    {
      wr[i] = hr[i+i*n];
      wi[i] = hi[i+i*n];
    }
  }

  en = igh;
  tr = 0.0;
  ti = 0.0;
  itn = 30 * n;
/*
  Search for next eigenvalue.
*/
  its = 0;
  enm1 = en - 1;
/*
  Look for single small sub-diagonal element.
*/
  while ( true )
  {
    for ( l = en; low <= l; l-- )
    {
      if ( l == low )
      {
        break;
      }

      tst1 = fabs ( hr[l-1+(l-1)*n] ) + fabs ( hi[l-1+(l-1)*n] ) + fabs ( hr[l+l*n] ) 
        + fabs ( hi[l+l*n] );
      tst2 = tst1 + fabs ( hr[l+(l-1)*n] ) + fabs ( hi[l+(l-1)*n] );

      if ( tst2 == tst1 )
      {
        break;
      }
    }
/*
  A root found.
*/
    if ( l == en )
    {
      hr[en+en*n] = hr[en+en*n] + tr;
      wr[en] = hr[en+en*n];
      hi[en+en*n] = hi[en+en*n] + ti;
      wi[en] = hi[en+en*n];
      en = enm1;
      if ( en < low )
      {
        break;
      }
      its = 0;
      enm1 = en - 1;
      continue;
    }

    if ( itn == 0 )
    {
      ierr = en;
      return ierr;
    }
/*
  Form shift.
*/
    if ( its == 10 || its == 20 )
    {
      sr = fabs ( hr[en+enm1*n] ) + fabs ( hr[enm1+(en-2)*n] );
      si = fabs ( hi[en+enm1*n] ) + fabs ( hi[enm1+(en-2)*n] );
    }
    else
    {
      sr = hr[en+en*n];
      si = hi[en+en*n];
      xr = hr[enm1+en*n] * hr[en+enm1*n] - hi[enm1+en*n] * hi[en+enm1*n];
      xi = hr[enm1+en*n] * hi[en+enm1*n] + hi[enm1+en*n] * hr[en+enm1*n];

      if ( xr != 0.0 || xi != 0.0 )
      {
        yr = ( hr[enm1+enm1*n] - sr ) / 2.0;
        yi = ( hi[enm1+enm1*n] - si ) / 2.0;
        csroot ( yr * yr - yi * yi + xr, 2.0 * yr * yi + xi, &zzr, &zzi );

        if ( yr * zzr + yi * zzi < 0.0 )
        {
          zzr = - zzr;
          zzi = - zzi;
        }
        cdiv ( xr, xi, yr + zzr, yi + zzi, &xr, &xi );
        sr = sr - xr;
        si = si - xi;
      }
    }

    for ( i = low; i <= en; i++ )
    {
      hr[i+i*n] = hr[i+i*n] - sr;
      hi[i+i*n] = hi[i+i*n] - si;
    }

    tr = tr + sr;
    ti = ti + si;
    its = its + 1;
    itn = itn - 1;
/*
  Look for two consecutive small sub-diagonal elements.
*/
    xr = fabs ( hr[enm1+enm1*n] ) + fabs ( hi[enm1+enm1*n] );
    yr = fabs ( hr[en+enm1*n] ) + fabs ( hi[en+enm1*n] );
    zzr = fabs ( hr[en+en*n] ) + fabs ( hi[en+en*n] );

    for ( m = enm1; l <= m; m-- )
    {
      if ( m == l )
      {
        break;
      }
      yi = yr;
      yr = fabs ( hr[m+(m-1)*n] ) + fabs ( hi[m+(m-1)*n] );
      xi = zzr;
      zzr = xr;
      xr = fabs ( hr[m-1+(m-1)*n] ) + fabs ( hi[m-1+(m-1)*n] );
      tst1 = zzr / yi * ( zzr + xr + xi );
      tst2 = tst1 + yr;
      if ( tst2 == tst1 )
      {
        break;
      }
    }
/*
  Triangular decomposition H=L*R.
*/
    for ( i = m + 1; i <= en; i++ )
    {
      xr = hr[i-1+(i-1)*n];
      xi = hi[i-1+(i-1)*n];
      yr = hr[i+(i-1)*n];
      yi = hi[i+(i-1)*n];
/*
  Interchange rows of HR and HI.
*/
      if ( fabs ( xr ) + fabs ( xi) < fabs ( yr ) + fabs ( yi ) )
      {
        for ( j = i - 1; j < n; j++ )
        {
          t           = hr[i-1+j*n];
          hr[i-1+j*n] = hr[i+j*n];
          hr[i+j*n]   = t;
          t           = hi[i-1+j*n];
          hi[i-1+j*n] = hi[i+j*n];
          hi[i+j*n]   = t;
        }

        cdiv ( xr, xi, yr, yi, &zzr, &zzi );
        wr[i] = 1.0;
      }
      else
      {
        cdiv ( yr, yi, xr, xi, &zzr, &zzi );
        wr[i] = - 1.0;
      }

      hr[i+(i-1)*n] = zzr;
      hi[i+(i-1)*n] = zzi;
      for ( j = i; j < n; j++ )
      {
        hr[i+j*n] = hr[i+j*n] - zzr * hr[i-1+j*n] + zzi * hi[i-1+j*n];
        hi[i+j*n] = hi[i+j*n] - zzr * hi[i-1+j*n] - zzi * hr[i-1+j*n];
      }
    }
/*
  Composition R*L=H.
*/
    for ( j = m + 1; j <= en; j++ )
    {
      xr = hr[j+(j-1)*n];
      xi = hi[j+(j-1)*n];
      hr[j+(j-1)*n] = 0.0;
      hi[j+(j-1)*n] = 0.0;
/*
  Interchange columns of HR, HI, ZR, and ZI.
*/
      if ( 0.0 < wr[j] )
      {
        for ( i = 0; i <= j; i++ )
        {
          t             = hr[i+(j-1)*n];
          hr[i+(j-1)*n] = hr[i+j*n];
          hr[i+j*n]     = t;
          t             = hi[i+(j-1)*n];
          hi[i+(j-1)*n] = hi[i+j*n];
          hi[i+j*n]     = t;
        }

        for ( i = low; i <= igh; i++ )
        {
          t             = zr[i+(j-1)*n];
          zr[i+(j-1)*n] = zr[i+j*n];
          zr[i+j*n]     = t;
          t             = zi[i+(j-1)*n];
          zi[i+(j-1)*n] = zi[i+j*n];
          zi[i+j*n]     = t;
        }
      }

      for ( i = 0; i <= j; i++ )
      {
        hr[i+(j-1)*n] = hr[i+(j-1)*n] + xr * hr[i+j*n] - xi * hi[i+j*n];
        hi[i+(j-1)*n] = hi[i+(j-1)*n] + xr * hi[i+j*n] + xi * hr[i+j*n];
      }
/*
  Accumulate transformations.
*/
      for ( i = low; i <= igh; i++ )
      {
        zr[i+(j-1)*n] = zr[i+(j-1)*n] + xr * zr[i+j*n] - xi * zi[i+j*n];
        zi[i+(j-1)*n] = zi[i+(j-1)*n] + xr * zi[i+j*n] + xi * zr[i+j*n];
      }
    }
  }
/*
  All roots found.
  Backsubstitute to find vectors of upper triangular form.
*/
  norm = 0.0;
  for ( i = 0; i < n; i++ )
  {
    for ( j = i; j < n; j++ )
    {
      tr = fabs ( hr[i+j*n] ) + fabs ( hi[i+j*n] );
      norm = r8_max ( norm, tr );
    }
  }

  hr[0+0*n] = norm;
  if ( n == 1 )
  {
    return ierr;
  }

  if ( norm == 0.0 )
  {
    return ierr;
  }

  for ( en = n - 1; 1 <= en; en-- )
  {
    xr = wr[en];
    xi = wi[en];
    hr[en+en*n] = 1.0;
    hi[en+en*n] = 0.0;
    enm1 = en - 1;

    for ( i = en - 1; 0 <= i; i-- )
    {
      zzr = 0.0;
      zzi = 0.0;
      for ( j = i + 1; j <= en; j++ )
      {
        zzr = zzr + hr[i+j*n] * hr[j+en*n] - hi[i+j*n] * hi[j+en*n];
        zzi = zzi + hr[i+j*n] * hi[j+en*n] + hi[i+j*n] * hr[j+en*n];
      }

      yr = xr - wr[i];
      yi = xi - wi[i];

      if ( yr == 0.0 && yi == 0.0 )
      {
        tst1 = norm;
        yr = tst1;

        while ( true )
        {
          yr = 0.01 * yr;
          tst2 = norm + yr;
          if ( tst2 <=  tst1 )
          {
            break;
          }
        }
      }
      cdiv ( zzr, zzi, yr, yi, &ar, &ai );
      hr[i+en*n] = ar;
      hi[i+en*n] = ai;
/*
  Overflow control.
*/
      tr = fabs ( hr[i+en*n] ) + fabs ( hi[i+en*n] );

      if ( tr != 0.0 )
      {
        tst1 = tr;
        tst2 = tst1 + 1.0 / tst1;

        if ( tst2 <= tst1 )
        {
          for ( j = i; j <= en; j++ )
          {
            hr[j+en*n] = hr[j+en*n] / tr;
            hi[j+en*n] = hi[j+en*n] / tr;
          }
        }
      }
    }
  }
/*
  End backsubstitution.
*/
  enm1 = n - 1;
/*
  Vectors of isolated roots.
*/
  for ( i = 0; i < n - 1; i++ )
  {
    if ( i < low || igh < i )
    {
      for ( j = i + 1; j < n; j++ )
      {
        zr[i+j*n] = hr[i+j*n];
        zi[i+j*n] = hi[i+j*n];
      }
    }
  }
/*
  Multiply by transformation matrix to give vectors of original full matrix.
*/
  for ( j = n - 1; low + 1 <= j; j-- )
  {
    m = i4_min ( j, igh );

    for ( i = low; i <= igh; i++ )
    {
      zzr = 0.0;
      zzi = 0.0;
      for ( k = low; k <= m; k++ )
      {
        zzr = zzr + zr[i+k*n] * hr[k+j*n] - zi[i+k*n] * hi[k+j*n];
        zzi = zzi + zr[i+k*n] * hi[k+j*n] + zi[i+k*n] * hr[k+j*n];
      }
      zr[i+j*n] = zzr;
      zi[i+j*n] = zzi;
    }
  }

  return ierr;
}
/******************************************************************************/

int comqr ( int n, int low, int igh, double hr[], double hi[], double wr[], 
  double wi[] )

/******************************************************************************/
/*
  Purpose:

    COMQR gets eigenvalues of a complex upper Hessenberg matrix.

  Discussion:

    COMQR finds the eigenvalues of a complex upper Hessenberg matrix by 
    the QR method.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing
    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.

    Input/output, double HR[N,N], HI[N,N].  On input, the real
    and imaginary parts of the complex upper Hessenberg matrix.  Their lower
    triangles below the subdiagonal contain information about the unitary
    transformations used in the reduction by CORTH, if performed.  On output,
    the upper Hessenberg portions of HR and HI have been destroyed.
    Therefore, they must be saved before calling COMQR if subsequent
    calculation of eigenvectors is to be performed.

    Output, double WR[N], WI[N], the real and imaginary parts of the
    eigenvalues.  If an error break; is made, the eigenvalues should be
    correct for indices IERR+1,...,N.

    Output, int COMQR, error flag.
    0, for normal return,
    J, if the limit of 30*N iterations is exhausted while the J-th
       eigenvalue is being sought.
*/
{
  double ai;
  double ar;
  int en;
  int enm1;
  int i;
  int ierr;
  int itn;
  int its;
  int j;
  int l;
  int ll;
  double norm;
  double si;
  double sr;
  double ti;
  double tr;
  double tst1;
  double tst2;
  double xi;
  double xr;
  double yi;
  double yr;
  double zzi;
  double zzr;

  ierr = 0;
/*
  Create real subdiagonal elements.
*/
  l = low + 1;

  for ( i = low + 1; i <= igh; i++ )
  {
    ll = i4_min ( i + 1, igh );

    if ( hi[i+(i-1)*n] != 0.0 )
    {
      norm = pythag ( hr[i+(i-1)*n], hi[i+(i-1)*n] );
      yr = hr[i+(i-1)*n] / norm;
      yi = hi[i+(i-1)*n] / norm;
      hr[i+(i-1)*n] = norm;
      hi[i+(i-1)*n] = 0.0;

      for ( j = i; j <= igh; j++ )
      {
        si = yr * hi[i+j*n] - yi * hr[i+j*n];
        hr[i+j*n] = yr * hr[i+j*n] + yi * hi[i+j*n];
        hi[i+j*n] = si;
      }

      for ( j = low; j <= ll; j++ )
      {
        si = yr * hi[j+i*n] + yi * hr[j+i*n];
        hr[j+i*n] = yr * hr[j+i*n] - yi * hi[j+i*n];
        hi[j+i*n] = si;
      }
    }
  }
/*
  Store roots isolated by CBAL.
*/
  for ( i = 0; i < n; i++ )
  {
    if ( i < low || igh < i )
    {
      wr[i] = hr[i+i*n];
      wi[i] = hi[i+i*n];
    }
  }

  en = igh;
  tr = 0.0;
  ti = 0.0;
  itn = 30 * n;
/*
  Search for next eigenvalue.
*/
  if ( en < low )
  {
    return ierr;
  }

  its = 0;
  enm1 = en - 1;
/*
  Look for single small sub-diagonal element.
*/
  while ( true )
  {
    for ( l = en; low <= l; l-- )
    {
      if ( l == low )
      {
        break;
      }
      tst1 = fabs ( hr[l-1+(l-1)*n] ) + fabs ( hi[l-1+(l-1)*n] ) + fabs ( hr[l+l*n] ) 
        + fabs ( hi[l+l*n] );
      tst2 = tst1 + fabs ( hr[l+(l-1)*n] );
      if ( tst2 == tst1 )
      {
        break;
      }
    } 
/*
  A root found.
*/
    if ( l == en )
    {
      wr[en] = hr[en+en*n] + tr;
      wi[en] = hi[en+en*n] + ti;
      en = enm1;
      if ( en < low )
      {
        break;
      }
      its = 0;
      enm1 = en - 1;
      continue;
    }

    if ( itn == 0 )
    {
      ierr = en;
      break;
    }

    if ( its == 10 || its == 20 )
    {
      sr = fabs ( hr[en+enm1*n] ) + fabs ( hr[enm1+(en-2)*n] );
      si = 0.0;
    }
    else
    {
      sr = hr[en+en*n];
      si = hi[en+en*n];
      xr = hr[enm1+en*n] * hr[en+enm1*n];
      xi = hi[enm1+en*n] * hr[en+enm1*n];

      if ( xr != 0.0 || xi != 0.0 )
      {
        yr = ( hr[enm1+enm1*n] - sr ) / 2.0;
        yi = ( hi[enm1+enm1*n] - si ) / 2.0;

        ar = yr * yr - yi * yi + xr;
        ai = 2.0 * yr * yi + xi;
        csroot ( ar, ai, &zzr, &zzi );

        if ( yr * zzr + yi * zzi < 0.0 )
        {
          zzr = - zzr;
          zzi = - zzi;
        }

        ar = yr + zzr;
        ai = yi + zzi;
        cdiv ( xr, xi, ar, ai, &xr, &xi );
        sr = sr - xr;
        si = si - xi;
      }
    }

    for ( i = low; i <= en; i++ )
    {
      hr[i+i*n] = hr[i+i*n] - sr;
      hi[i+i*n] = hi[i+i*n] - si;
    }

    tr = tr + sr;
    ti = ti + si;
    its = its + 1;
    itn = itn - 1;
/*
  Reduce to triangle (rows).
*/
    for ( i = l + 1; i <= en; i++ )
    {
      sr = hr[i+(i-1)*n];
      hr[i+(i-1)*n] = 0.0;
      norm = pythag ( pythag ( hr[i-1+(i-1)*n], hi[i-1+(i-1)*n] ), sr );
      xr = hr[i-1+(i-1)*n] / norm;
      wr[i-1] = xr;
      xi = hi[i-1+(i-1)*n] / norm;
      wi[i-1] = xi;
      hr[i-1+(i-1)*n] = norm;
      hi[i-1+(i-1)*n] = 0.0;
      hi[i+(i-1)*n] = sr / norm;

      for ( j = i; j <= en; j++ )
      {
        yr = hr[i-1+j*n];
        yi = hi[i-1+j*n];
        zzr = hr[i+j*n];
        zzi = hi[i+j*n];
        hr[i-1+j*n] = xr * yr + xi * yi + hi[i+(i-1)*n] * zzr;
        hi[i-1+j*n] = xr * yi - xi * yr + hi[i+(i-1)*n] * zzi;
        hr[i+j*n] = xr * zzr - xi * zzi - hi[i+(i-1)*n] * yr;
        hi[i+j*n] = xr * zzi + xi * zzr - hi[i+(i-1)*n] * yi;
      }
    }

    si = hi[en+en*n];

    if ( si != 0.0 )
    {
      norm = pythag ( hr[en+en*n], si );
      sr = hr[en+en*n] / norm;
      si = si / norm;
      hr[en+en*n] = norm;
      hi[en+en*n] = 0.0;
    }
/*
  Inverse operation (columns).
*/
    for ( j = l + 1; j <= en; j++ )
    {
      xr = wr[j-1];
      xi = wi[j-1];

      for ( i = l; i <= j; i++ )
      {
        yr = hr[i+(j-1)*n];
        yi = 0.0;
        zzr = hr[i+j*n];
        zzi = hi[i+j*n];
        if ( i != j )
        {
          yi = hi[i+(j-1)*n];
          hi[i+(j-1)*n] = xr * yi + xi * yr + hi[j+(j-1)*n] * zzi;
        }
        hr[i+(j-1)*n] = xr * yr - xi * yi + hi[j+(j-1)*n] * zzr;
        hr[i+j*n] = xr * zzr + xi * zzi - hi[j+(j-1)*n] * yr;
        hi[i+j*n] = xr * zzi - xi * zzr - hi[j+(j-1)*n] * yi;
      }
    }

    if ( si != 0.0 )
    {
      for ( i = l; i <= en; i++ )
      {
        yr = hr[i+en*n];
        yi = hi[i+en*n];
        hr[i+en*n] = sr * yr - si * yi;
        hi[i+en*n] = sr * yi + si * yr;
      }
    }
  }

  return ierr;
}
/******************************************************************************/

int comqr2 ( int n, int low, int igh, double ortr[], double orti[], 
  double hr[], double hi[], double wr[], double wi[], double zr[], 
  double zi[] )

/******************************************************************************/
/*
  Purpose:

    COMQR2 gets eigenvalues/vectors of a complex upper Hessenberg matrix.

  Discussion:

    COMQR2 finds the eigenvalues and eigenvectors
    of a complex upper Hessenberg matrix by the QR
    method.  The eigenvectors of a complex general matrix
    can also be found if CORTH has been used to reduce
    this general matrix to Hessenberg form.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing
    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.

    Input/output, double ORTR[N], ORTI[N].  On input, information 
    about the unitary transformations used in the reduction by CORTH, if 
    performed.  If the eigenvectors of the Hessenberg matrix are desired, set 
    ORTR(J) and ORTI(J) to 0.0 for these elements.  On output, these arrays
    have been overwritten.

    Input/output, double HR[N,N], HI[N,N].  On input, the real and 
    imaginary parts of the complex upper Hessenberg matrix.  Their lower 
    triangles below the subdiagonal contain further information about the
    transformations which were used in the reduction by CORTH, if performed.
    If the eigenvectors of the Hessenberg matrix are desired, these elements
    may be arbitrary.

    Output, double WR[N], WI[N], the real and imaginary parts of the
    eigenvalues.  If an error exit is made, the eigenvalues should be
    correct for indices IERR+1,...,N.

    Output, double ZR[N,N], ZI[N,N], the real and imaginary parts of
    the eigenvectors.  The eigenvectors are unnormalized.  If an error exit
    is made, none of the eigenvectors has been found.

    Output, int COMQR2, error flag.
    0, for normal return,
    J, if the limit of 30*N iterations is exhausted while the J-th
      eigenvalue is being sought.
*/
{
  double ai;
  double ar;
  int en;
  int enm1;
  int i;
  int ierr;
  int itn;
  int its;
  int j;
  int k;
  int l;
  int ll;
  int m;
  int nn;
  double norm;
  double si;
  double sr;
  double ti;
  double tr;
  double tst1;
  double tst2;
  double xi;
  double xr;
  double yi;
  double yr;
  double zzi;
  double zzr;

  ierr = 0;
/*
  Initialize eigenvector matrix.
*/
  r8mat_identity ( n, zr );
  r8mat_zeros ( n, n, zi );
/*
  Form the matrix of accumulated transformations from the information
  left by CORTH.
*/
  if ( 0 < igh - low - 1 )
  {
    for ( i = igh - 1; low + 1 <= i; i-- )
    {
      if ( ortr[i] == 0.0 && orti[i] == 0.0 )
      {
        continue;
      }

      if ( hr[i+(i-1)*n] == 0.0 && hi[i+(i-1)*n] == 0.0 )
      {
        continue;
      }
/*
  Norm below is negative of H formed in CORTH.
*/
      norm = hr[i+(i-1)*n] * ortr[i] + hi[i+(i-1)*n] * orti[i];

      for ( k = i + 1; k <= igh; k++ )
      {
        ortr[k] = hr[k+(i-1)*n];
        orti[k] = hi[k+(i-1)*n];
      }

      for ( j = i; j <= igh; j++ )
      {
        sr = 0.0;
        si = 0.0;
        for ( k = i; k <= igh; k++ )
        {
          sr = sr + ortr[k] * zr[k+j*n] + orti[k] * zi[k+j*n];
          si = si + ortr[k] * zi[k+j*n] - orti[k] * zr[k+j*n];
        }
        sr = sr / norm;
        si = si / norm;

        for ( k = i; k <= igh; k++ )
        {
          zr[k+j*n] = zr[k+j*n] + sr * ortr[k] - si * orti[k];
          zi[k+j*n] = zi[k+j*n] + sr * orti[k] + si * ortr[k];
        }
      }
    }
  }
/*
  Create real subdiagonal elements.
*/
  if ( 0 <= igh - low - 1 )
  {
    l = low + 1;

    for ( i = low + 1; i <= igh; i++ )
    {
      ll = i4_min ( i + 1, igh );

      if ( hi[i+(i-1)*n] == 0.0 )
      {
        continue;
      }

      norm = pythag ( hr[i+(i-1)*n], hi[i+(i-1)*n] );
      yr = hr[i+(i-1)*n] / norm;
      yi = hi[i+(i-1)*n] / norm;
      hr[i+(i-1)*n] = norm;
      hi[i+(i-1)*n] = 0.0;

      for ( j = i; j < n; j++ )
      {
        si =        yr * hi[i+j*n] - yi * hr[i+j*n];
        hr[i+j*n] = yr * hr[i+j*n] + yi * hi[i+j*n];
        hi[i+j*n] = si;
      }

      for ( j = 0; j <= ll; j++ )
      {
        si =        yr * hi[j+i*n] + yi * hr[j+i*n];
        hr[j+i*n] = yr * hr[j+i*n] - yi * hi[j+i*n];
        hi[j+i*n] = si;
      }

      for ( j = low; j <= igh; j++ )
      {
        si =        yr * zi[j+i*n] + yi * zr[j+i*n];
        zr[j+i*n] = yr * zr[j+i*n] - yi * zi[j+i*n];
        zi[j+i*n] = si;
      }
    }
  }
/*
  Store roots isolated by CBAL.
*/
  for ( i = 0; i < n; i++ )
  {
    if ( i < low || igh < i )
    {
      wr[i] = hr[i+i*n];
      wi[i] = hi[i+i*n];
    }
  }

  en = igh;
  tr = 0.0;
  ti = 0.0;
  itn = 30 * n;
/*
  Search for next eigenvalue.
*/
  its = 0;
  enm1 = en - 1;
/*
  Look for single small sub-diagonal element.
*/
  while ( true )
  {
    for ( l = en; low + 1 <= l; l-- )
    {
      tst1 = fabs ( hr[l-1+(l-1)*n] ) + fabs ( hi[l-1+(l-1)*n] ) 
        + fabs ( hr[l+l*n] ) + fabs ( hi[l+l*n] );
      tst2 = tst1 + fabs ( hr[l+(l-1)*n] );
      if ( tst2 == tst1 )
      {
        break;
      }
    }
/*
  A root found.
*/
    if ( l == en )
    {
      hr[en+en*n] = hr[en+en*n] + tr;
      wr[en] = hr[en+en*n];
      hi[en+en*n] = hi[en+en*n] + ti;
      wi[en] = hi[en+en*n];
      en = enm1;
      if ( en < low )
      {
        break;
      }
      its = 0;
      enm1 = en - 1;
      continue;
    }

    if ( itn == 0 )
    {
      ierr = en;
      return ierr;
    }

    if ( its == 10 || its == 20 )
    {
      sr = fabs ( hr[en+enm1*n] ) + fabs ( hr[enm1+(en-2)*n] );
      si = 0.0;
    }
    else
    {
      sr = hr[en+en*n];
      si = hi[en+en*n];
      xr = hr[enm1+en*n] * hr[en+enm1*n];
      xi = hi[enm1+en*n] * hr[en+enm1*n];

      if ( xr != 0.0 || xi != 0.0 )
      {
        yr = ( hr[enm1+enm1*n] - sr ) / 2.0;
        yi = ( hi[enm1+enm1*n] - si ) / 2.0;

        csroot ( yr * yr - yi * yi + xr, 2.0 * yr * yi + xi, &zzr, &zzi );

        if ( yr * zzr + yi * zzi < 0.0 )
        {
          zzr = - zzr;
          zzi = - zzi;
        }

        cdiv ( xr, xi, yr + zzr, yi + zzi, &xr, &xi );
        sr = sr - xr;
        si = si - xi;
      }
    }

    for ( i = low; i <= en; i++ )
    {
      hr[i+i*n] = hr[i+i*n] - sr;
      hi[i+i*n] = hi[i+i*n] - si;
    }

    tr = tr + sr;
    ti = ti + si;
    its = its + 1;
    itn = itn - 1;
/*
  Reduce to triangle (rows).
*/
    for ( i = l + 1; i <= en; i++ )
    {
      sr = hr[i+(i-1)*n];
      hr[i+(i-1)*n] = 0.0;
      norm = pythag ( pythag ( hr[i-1+(i-1)*n], hi[i-1+(i-1)*n] ), sr );
      xr = hr[i-1+(i-1)*n] / norm;
      wr[i-1] = xr;
      xi = hi[i-1+(i-1)*n] / norm;
      wi[i-1] = xi;
      hr[i-1+(i-1)*n] = norm;
      hi[i-1+(i-1)*n] = 0.0;
      hi[i+(i-1)*n] = sr / norm;

      for ( j = i; j < n; j++ )
      {
        yr = hr[i-1+j*n];
        yi = hi[i-1+j*n];
        zzr = hr[i+j*n];
        zzi = hi[i+j*n];
        hr[i-1+j*n] = xr * yr + xi * yi + hi[i+(i-1)*n] * zzr;
        hi[i-1+j*n] = xr * yi - xi * yr + hi[i+(i-1)*n] * zzi;
        hr[i+j*n] = xr * zzr - xi * zzi - hi[i+(i-1)*n] * yr;
        hi[i+j*n] = xr * zzi + xi * zzr - hi[i+(i-1)*n] * yi;
      }
    }

    si = hi[en+en*n];

    if ( si != 0.0 )
    {
      norm = pythag ( hr[en+en*n], si );
      sr = hr[en+en*n] / norm;
      si = si / norm;
      hr[en+en*n] = norm;
      hi[en+en*n] = 0.0;

      for ( j = en + 1; j < n; j++ )
      {
        yr = hr[en+j*n];
        yi = hi[en+j*n];
        hr[en+j*n] = sr * yr + si * yi;
        hi[en+j*n] = sr * yi - si * yr;
      }
    }
/*
  Inverse operation (columns).
*/
    for ( j = l + 1; j <= en; j++ )
    {
      xr = wr[j-1];
      xi = wi[j-1];

      for ( i = 0; i <= j; i++ )
      {
        yr = hr[i+(j-1)*n];
        yi = 0.0;
        zzr = hr[i+j*n];
        zzi = hi[i+j*n];

        if ( i != j )
        {
          yi = hi[i+(j-1)*n];
          hi[i+(j-1)*n] = xr * yi + xi * yr + hi[j+(j-1)*n] * zzi;
        }
        hr[i+(j-1)*n] = xr * yr - xi * yi + hi[j+(j-1)*n] * zzr;
        hr[i+j*n] = xr * zzr + xi * zzi - hi[j+(j-1)*n] * yr;
        hi[i+j*n] = xr * zzi - xi * zzr - hi[j+(j-1)*n] * yi;
      }

      for ( i = low; i <= igh; i++ )
      {
        yr = zr[i+(j-1)*n];
        yi = zi[i+(j-1)*n];
        zzr = zr[i+j*n];
        zzi = zi[i+j*n];
        zr[i+(j-1)*n] = xr * yr - xi * yi + hi[j+(j-1)*n] * zzr;
        zi[i+(j-1)*n] = xr * yi + xi * yr + hi[j+(j-1)*n] * zzi;
        zr[i+j*n] = xr * zzr + xi * zzi - hi[j+(j-1)*n] * yr;
        zi[i+j*n] = xr * zzi - xi * zzr - hi[j+(j-1)*n] * yi;
      }
    }

    if ( si != 0.0 )
    {
      for ( i = 0; i <= en; i++ )
      {
        yr = hr[i+en*n];
        yi = hi[i+en*n];
        hr[i+en*n] = sr * yr - si * yi;
        hi[i+en*n] = sr * yi + si * yr;
      }

      for ( i = low; i <= igh; i++ )
      {
        yr = zr[i+en*n];
        yi = zi[i+en*n];
        zr[i+en*n] = sr * yr - si * yi;
        zi[i+en*n] = sr * yi + si * yr;
      }
    }
  }
/*
  All roots found.
  Backsubstitute to find vectors of upper triangular form.
*/
  norm = 0.0;
  for ( i = 0; i < n; i++ )
  {
    for ( j = i; j < n; j++ )
    {
      tr = fabs ( hr[i+j*n] ) + fabs ( hi[i+j*n] );
      norm = r8_max ( norm, tr );
    }
  }

  if ( n == 1 )
  {
    return ierr;
  }

  if ( norm == 0.0 )
  {
    return ierr;
  }

  for ( nn = 1; nn < n; nn++ )
  {
    en = n + 2 - nn;
    xr = wr[en];
    xi = wi[en];
    hr[en+en*n] = 1.0;
    hi[en+en*n] = 0.0;
    enm1 = en - 1;

    for ( i = en - 1; 0 <= i; i-- )
    {
      zzr = 0.0;
      zzi = 0.0;
      for ( j = i + 1; j <= en; j++ )
      {
        zzr = zzr + hr[i+j*n] * hr[j+en*n] - hi[i+j*n] * hi[j+en*n];
        zzi = zzi + hr[i+j*n] * hi[j+en*n] + hi[i+j*n] * hr[j+en*n];
      }
      yr = xr - wr[i];
      yi = xi - wi[i];

      if ( yr == 0.0 && yi == 0.0 )
      {
        tst1 = norm;
        yr = tst1;

        while ( true )
        {
          yr = 0.01 * yr;
          tst2 = norm + yr;
          if ( tst2 <= tst1 )
          {
            break;
          }
        }
      }

      cdiv ( zzr, zzi, yr, yi, &ar, &ai );
      hr[i+en*n] = ar;
      hi[i+en*n] = ai;
/*
  Overflow control.
*/
      tr = fabs ( hr[i+en*n] ) + fabs ( hi[i+en*n] );

      if ( tr != 0.0 )
      {
        tst1 = tr;
        tst2 = tst1 + 1.0 / tst1;

        if ( tst2 <= tst1 )
        {
          for ( j = i; j <= en; j++ )
          {
            hr[j+en*n] = hr[j+en*n] / tr;
            hi[j+en*n] = hi[j+en*n] / tr;
          }
        }
      }
    }
  }
/*
  End backsubstitution.
*/
  enm1 = n - 1;
/*
  Vectors of isolated roots.
*/
  for ( i = 0; i < n - 1; i++ )
  {
    if ( i < low || igh < i )
    {
      for ( j = i + 1; j < n; j++ )
      {
        zr[i+j*n] = hr[i+j*n];
        zi[i+j*n] = hi[i+j*n];
      }
    }
  }
/*
  Multiply by transformation matrix to give vectors of original full matrix.
*/
  for ( j = n - 1; low + 1 <= j; j-- )
  {
    m = i4_min ( j, igh );

    for ( i = low; i <= igh; i++ )
    {
      zzr = 0.0;
      zzi = 0.0;
      for ( k = low; k <= m; k++ )
      {
        zzr = zzr + zr[i+k*n] * hr[k+j*n] - zi[i+k*n] * hi[k+j*n];
        zzi = zzi + zr[i+k*n] * hi[k+j*n] + zi[i+k*n] * hr[k+j*n];
      }
      zr[i+j*n] = zzr;
      zi[i+j*n] = zzi;
    }
  }

  return ierr;
}

/****=**************************************************************************/

void cortb ( int n, int low, int igh, double ar[], double ai[], double ortr[],
  double orti[], int m, double zr[], double zi[] )

/******************************************************************************/
/*
  Purpose:

    CORTB determines eigenvectors by undoing the CORTH transformation.

  Discussion:

    CORTB forms the eigenvectors of a complex general
    matrix by back transforming those of the corresponding
    upper Hessenberg matrix determined by CORTH.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing
    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH to the order
    of the matrix.

    Input, double AR[N,IGH], AI[N,IGH], information about the 
    unitary transformations used in the reduction by CORTH in their strict 
    lower triangles.

    Input/output, double ORTR[IGH], ORTI[IGH].  On input, further 
    information about the transformations used in the reduction by CORTH.  On 
    output, ORTR and ORTI have been further altered.

    Input, int M, the number of columns of ZR and ZI to be 
    back transformed.

    Input/output, double ZR[N,M], ZI[N,M].  On input, the real and 
    imaginary parts of the eigenvectors to be back transformed.  On output, 
    the real and imaginary parts of the transformed eigenvectors.
*/
{
  double gi;
  double gr;
  double h;
  int i;
  double ii = 0.0 ;
  double ir;
  int j;
  int mp;
  double ri;
  double rr;

  if ( m == 0 )
  {
    return;
  }

  if ( igh - 1 < low + 1 )
  {
    return;
  }

  for ( mp = igh - 1; low + 1 <= mp; mp-- )
  {
    if ( ar[mp+(mp-1)*n] != 0.0 || ai[mp+(mp-1)*n] != 0.0 )
    {
      h = ar[mp+(mp-1)*n] * ortr[mp] + ai[mp+(mp-1)*n] * orti[mp];
      for ( i = mp + 1; i <= igh; i++ )
      {
        ortr[i] = ar[i+(mp-1)*n];
        orti[i] = ai[i+(mp-1)*n];
      }
      for ( j = 0; j < m; j++ )
      {
        rr = 0.0;
        ri = 0.0;
        for ( i = mp; i <= igh; i++ )
        {
          rr = rr + ortr[i] * zr[i+j*n];
          ii = ii + orti[i] * zi[i+j*n];
        }
        gr = ( rr + ii ) / h;

        ri = 0.0;
        ir = 0.0;
        for ( i = mp; i <= igh; i++ )
        {
          ri = ri + ortr[i] * zi[i+j*n];
          ir = ir + orti[i] * zr[i+j*n];
        }
        gi = ( ri - ir ) / h;

        for ( i = mp; i <= igh; i++ )
        {
          zr[i+j*n] = zr[i+j*n] + gr * ortr[i] - gi * orti[i];
          zi[i+j*n] = zi[i+j*n] + gr * orti[i] + gi * ortr[i];
        }
      }
    }
  }

  return;
}
/******************************************************************************/

void corth ( int n, int low, int igh, double ar[], double ai[], double ortr[], 
  double orti[] )

/******************************************************************************/
/*
  Purpose:

    CORTH transforms a complex general matrix to upper Hessenberg form.

  Discussion:

    CORTH is given a complex general matrix and reduces a submatrix situated 
    in rows and columns LOW through IGH to upper Hessenberg form by
    unitary similarity transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing
    routine CBAL.  If CBAL is not used, set LOW = 1 and IGH = N.

    Input/output, double AR[N,N], AI[N,N].  On input, the real and
    imaginary parts of the complex input matrix.  On output, the real and 
    imaginary parts of the Hessenberg matrix.  Information about the unitary
    transformations used in the reduction is stored in the remaining
    triangles under the Hessenberg matrix.

    Output, double ORTR[IGH], ORTI[IGH], further information about 
    the transformations.
*/
{
  double f;
  double fi;
  double fr;
  double g;
  double h;
  int i;
  int j;
  int m;
  double scale;

  if ( igh - 1 < low + 1 )
  {
    return;
  }

  for ( m = low + 1; m <= igh - 1; m++ )
  {
    h = 0.0;
    ortr[m] = 0.0;
    orti[m] = 0.0;
    scale = 0.0;
/*
  Scale column.
*/
    for ( i = m; i <= igh; i++ )
    {
      scale = scale + fabs ( ar[i+(m-1)*n] ) + fabs ( ai[i+(m-1)*n] );
    }

    if ( scale == 0.0 )
    {
      continue;
    }

    for ( i = igh; m <= i; i-- )
    {
      ortr[i] = ar[i+(m-1)*n] / scale;
      orti[i] = ai[i+(m-1)*n] / scale;
      h = h + ortr[i] * ortr[i] + orti[i] * orti[i];
    }

    g = sqrt ( h );
    f = pythag ( ortr[m], orti[m] );

    if ( f != 0.0 )
    {
      h = h + f * g;
      g = g / f;
      ortr[m] = ( 1.0 + g ) * ortr[m];
      orti[m] = ( 1.0 + g ) * orti[m];
    }
    else
    {
      ortr[m] = g;
      ar[m+(m-1)*n] = scale;
    }
/*
  Form (I-(U*Ut)/h) * A.
*/
    for ( j = m; j < n; j++ )
    {
      fr = 0.0;
      fi = 0.0;
      for ( i = igh; m <= i; i-- )
      {
        fr = fr + ortr[i] * ar[i+j*n] + orti[i] * ai[i+j*n];
        fi = fi + ortr[i] * ai[i+j*n] - orti[i] * ar[i+j*n];
      }
      fr = fr / h;
      fi = fi / h;

      for ( i = m; i <= igh; i++ )
      {
        ar[i+j*n] = ar[i+j*n] - fr * ortr[i] + fi * orti[i];
        ai[i+j*n] = ai[i+j*n] - fr * orti[i] - fi * ortr[i];
      }
    }
/*
  Form (I-(U*Ut)/h) * A * (I-(U*Ut)/h)
*/
    for ( i = 0; i <= igh; i++ )
    {
      fr = 0.0;
      fi = 0.0;
      for ( j = igh; m <= j; j-- )
      {
        fr = fr + ortr[j] * ar[i+j*n] - orti[j] * ai[i+j*n];
        fi = fi + ortr[j] * ai[i+j*n] + orti[j] * ar[i+j*n];
      }
      fr = fr / h;
      fi = fi / h;

      for ( j = m; j <= igh; j++ )
      {
        ar[i+j*n] = ar[i+j*n] - fr * ortr[j] - fi * orti[j];
        ai[i+j*n] = ai[i+j*n] + fr * orti[j] - fi * ortr[j];
      }
    }

    ortr[m] = scale * ortr[m];
    orti[m] = scale * orti[m];
    ar[m+(m-1)*n] = - g * ar[m+(m-1)*n];
    ai[m+(m-1)*n] = - g * ai[m+(m-1)*n];
  }

  return;
}
/******************************************************************************/

void csroot ( double xr, double xi, double *yr, double *yi )

/******************************************************************************/
/*
  Purpose:

    CSROOT computes the complex square root of a complex quantity.

  Discussion:

    The branch of the square function is chosen so that
      0.0 <= YR
    and
      sign ( YI ) == sign ( XI )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, double XR, XI, the real and imaginary parts of the 
    quantity whose square root is desired.

    Output, double *YR, *YI, the real and imaginary parts of the 
    square root.
*/
{
  double s;
  double ti;
  double tr;

  tr = xr;
  ti = xi;
  s = sqrt ( 0.5 * ( pythag ( tr, ti ) + fabs ( tr ) ) );

  if ( 0.0 <= tr )
  {
    *yr = s;
  }

  if ( ti < 0.0 )
  {
    s = -s;
  }

  if ( tr <= 0.0 )
  {
    *yi = s;
  }

  if ( tr < 0.0 )
  {
    *yr = 0.5 * ( ti / *yi );
  }
  else if ( 0.0 < tr )
  {
    *yi = 0.5 * ( ti / *yr );
  }
  return;
}
/******************************************************************************/

void elmbak ( int n, int low, int igh, double a[], int ind[], int m, double z[] )

/******************************************************************************/
/*
  Purpose:

    ELMBAK determines eigenvectors by undoing the ELMHES transformation.

  Discussion:

    ELMBAK forms the eigenvectors of a real general
    matrix by back transforming those of the corresponding
    upper Hessenberg matrix determined by ELMHES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, integers determined by the balancing
    routine BALANC.  If BALANC has not been used, set LOW = 1 and
    IGH equal to the order of the matrix.

    Input, double A[N,IGH], the multipliers which were used in the
    reduction by ELMHES in its lower triangle below the subdiagonal.

    Input, int IND[IGH]), information on the rows and columns
    interchanged in the reduction by ELMHES.

    Input, int M, the number of columns of Z to be back
    transformed.

    Input/output, double Z[N,M].  On input, the real and imaginary 
    parts of the eigenvectors to be back transformed.  On output, the real and
    imaginary parts of the transformed eigenvectors.
*/
{
  int i;
  int j;
  int la;
  int mm;
  int mp;
  double t;
  double x;

  if ( m == 0 )
  {
    return;
  }

  la = igh - 1;

  if ( la < low + 1 )
  {
    return;
  }

  for ( mm = low + 1; mm <= la; mm++ )
  {
    mp = low + igh - mm;

    for ( i = mp + 1; i < igh; i++ )
    {
      x = a[i+(mp-1)*n];
      if ( x != 0.0 )
      {
        for ( j = 0; j < m; j++ )
        {
          z[i+j*n] = z[i+j*n] + x * z[mp+j*n];
        }
      }
    }

    i = ind[mp];

    if ( i != mp )
    {
      for ( j = 0; j < m; j++ )
      {
        t         = z[i+j*n];
        z[i+j*n]  = z[mp+j*n];
        z[mp+j*n] = t;
      }
    }
  }

  return;
}
/******************************************************************************/

void elmhes ( int n, int low, int igh, double a[], int ind[] )

/******************************************************************************/
/*
  Purpose:

    ELMHES transforms a real general matrix to upper Hessenberg form.

  Discussion:

    ELMHES is given a real general matrix and reduces a submatrix
    situated in rows and columns LOW through IGH to upper Hessenberg
    form by stabilized elementary similarity transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Martin, James Wilkinson,
    ELMHES,
    Numerische Mathematik,
    Volume 12, pages 349-368, 1968.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing 
    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.

    Input/output, double A[N,N].  On input, the matrix to be 
    reduced.  On output, the Hessenberg matrix.  The multipliers
    which were used in the reduction are stored in the
    remaining triangle under the Hessenberg matrix.

    Output, int IND[N], contains information on the rows
    and columns interchanged in the reduction.  Only elements LOW through
    IGH are used.
*/
{
  int i;
  int j;
  int m;
  double t;
  double x;
  double y;

  for ( m = low + 1; m <= igh - 1; m++ )
  {
    x = 0.0;
    i = m;

    for ( j = m; j <= igh; j++ )
    {
      if ( fabs ( x ) < fabs ( a[j+(m-1)*n] ) )
      {
        x = a[j+(m-1)*n];
        i = j;
      }
    }

    ind[m] = i;
/*
  Interchange rows and columns of the matrix.
*/
    if ( i != m )
    {
      for ( j = m - 1; j < n; j++ )
      {
        t        = a[i+j*n];
        a[i+j*n] = a[m+j*n];
        a[m+j*n] = t;
      }

      for ( j = 0; j <= igh; j++ )
      {
        t        = a[j+i*n];
        a[j+i*n] = a[j+m*n];
        a[j+m*n] = t;
      }
    }

    if ( x != 0.0 )
    {
      for ( i = m + 1; i <= igh; i++ )
      {
        y = a[i+(m-1)*n];

        if ( y != 0.0 )
        {
          y = y / x;
          a[i+(m-1)*n] = y;

          for ( j = m; j < n; j++ )
          {
            a[i+j*n] = a[i+j*n] - y * a[m+j*n];
          }
          for ( j = 0; j <= igh; j++ )
          {
            a[j+m*n] = a[j+m*n] + y * a[j+i*n];
          }
        }
      }
    }
  }

  return;
}
/******************************************************************************/

void eltran ( int n, int low, int igh, double a[], int ind[], double z[] )

/******************************************************************************/
/*
  Purpose:

    ELTRAN accumulates similarity transformations used by ELMHES.

  Discussion:

    ELTRAN accumulates the stabilized elementary similarity transformations 
    used in the reduction of a real general matrix to upper Hessenberg form 
    by ELMHES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Peters, James Wilkinson,
    ELMTRANS,
    Numerische Mathematik,
    Volume 16, pages 181-204, 1970.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing 
    routine BALANC.  If BALANC has not been used, set LOW = 0, IGH = N - 1.

    Input, double A[N,IGH], the multipliers which were used in the
    reduction by ELMHES in its lower triangle below the subdiagonal.

    Input, int IND[IGH], information on the rows and columns
    interchanged in the reduction by ELMHES.

    Output, double Z[N,N], the transformation matrix produced in the
    reduction by ELMHES.
*/
{
  int i;
  int j;
  int mp;
/*
  Initialize Z to the identity matrix.
*/
  r8mat_identity ( n, z );

  if ( igh < low + 2 )
  {
    return;
  }

  for ( mp = igh - 1; low + 1 <= mp; mp-- )
  {
    for ( i = mp + 1; i <= igh; i++ )
    {
      z[i+mp*n] = a[i+(mp-1)*n];
    }
    i = ind[mp];

    if ( i != mp )
    {
      for ( j = mp; j <= igh; j++ )
      {
        z[mp+j*n] = z[i+j*n];
      }
      z[i+mp*n] = 1.0;
      for ( j = mp + 1; j <= igh; j++ )
      {
        z[i+j*n] = 0.0;
      }
    }
  }

  return;
}
/******************************************************************************/

int figi ( int n, double t[], double d[], double e[], double e2[] )

/******************************************************************************/
/*
  Purpose:

    FIGI transforms a real nonsymmetric tridiagonal matrix to symmetric form.

  Discussion:

    FIGI is given a nonsymmetric tridiagonal matrix such that the products
    of corresponding pairs of off-diagonal elements are all
    non-negative.  It reduces the matrix to a symmetric
    tridiagonal matrix with the same eigenvalues.  If, further,
    a zero product only occurs when both factors are zero,
    the reduced matrix is similar to the original matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double T[N,3] contains the input matrix.  Its subdiagonal
    is stored in the last N-1 positions of the first column, its diagonal in
    the N positions of the second column, and its superdiagonal in the
    first N-1 positions of the third column.  T(1,1) and T(N,3) are arbitrary.

    Output, double D[N], the diagonal elements of the symmetric 
    matrix.

    Output, double E[N], contains the subdiagonal elements of
    the symmetric matrix in E(2:N).  E(1) is not set.

    Output, double E2[N], the squares of the corresponding elements 
    of E.  E2 may coincide with E if the squares are not needed.

    Output, int FIGI, error flag.
    0, for normal return,
    N+I, if T(I,1) * T(I-1,3) is negative,
    -(3*N+I), if T(I,1) * T(I-1,3) is zero with one factor non-zero.  In
      this case, the eigenvectors of the symmetric matrix are not simply
      related to those of T and should not be sought.
*/
{
  int i;
  int ierr;

  ierr = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( 0 < i )
    {
      e2[i] = t[i+0*n] * t[i-1+2*n];

      if ( e2[i] < 0.0 )
      {
        ierr = n + i + 1;
        return ierr;
      }
      else if ( e2[i] == 0.0 )
      {
        if ( t[i+0*n] != 0.0 || t[i-1+2*n] != 0.0 )
        {
          ierr = - 3 * n - i - 1;
          return ierr;
        }
        e[i] = 0.0;
      }
      else
      {
        e[i] = sqrt ( e2[i] );
      }
    }
    else
    {
      e[i] = 0.0;
      e2[i] = 0.0;
    }
    d[i] = t[i+1*n];
  }

  return ierr;
}
/******************************************************************************/

int figi2 ( int n, double t[], double d[], double e[], double z[] )

/******************************************************************************/
/*
  Purpose:

    FIGI2 transforms a real nonsymmetric tridiagonal matrix to symmetric form.

  Discussion:

    FIGI2 is given a nonsymmetric tridiagonal matrix such that the products
    of corresponding pairs of off-diagonal elements are all
    non-negative, and zero only when both factors are zero.

    It reduces the matrix to a symmetric tridiagonal matrix
    using and accumulating diagonal similarity transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double T[N,3] contains the input matrix.  Its subdiagonal
    is stored in the last N-1 positions of the first column, its diagonal in
    the N positions of the second column, and its superdiagonal in the
    first N-1 positions of the third column.  T(1,1) and T(N,3) are arbitrary.

    Output, double D[N], the diagonal elements of the symmetric 
    matrix.

    Output, double E[N], contains the subdiagonal elements of the 
    symmetric matrix in E(2:N).  E(1) is not set.

    Output, double Z[N,N], contains the transformation matrix
    produced in the reduction.

    Output, int FIGI2, error flag.
    0, for normal return,
    N+I, if T(I,1) * T(I-1,3) is negative,
    2*N+I, if T(I,1) * T(I-1,3) is zero with one factor non-zero.
*/
{
  double h;
  int i;
  int ierr;

  ierr = 0;
/*
  Initialize Z to the identity matrix.
*/
  r8mat_identity ( n, z );

  for ( i = 0; i < n; i++ )
  {
    if ( i == 0 )
    {
      e[i] = 0.0;
    }
    else
    {
      h = t[i+0*n] * t[i-1+2*n];

      if ( h < 0.0 )
      {
        ierr = n + i + 1;
        return ierr;
      }
      else if ( h == 0 )
      {
        if ( t[i+0*n] != 0.0 || t[i-1+2*n] != 0.0 )
        {
          ierr = 2 * n + i + 1;
          return ierr;
        }

        e[i] = 0.0;
        z[i+i*n] = 1.0;
      }
      else
      {
        e[i] = sqrt ( h );
        z[i+i*n] = z[i-1+(i-1)*n] * e[i] / t[i-1+2*n];
      }
    }
    d[i] = t[i+1*n];
  }

  return ierr;
}
/******************************************************************************/

int hqr ( int n, int low, int igh, double h[], double wr[], double wi[] )

/******************************************************************************/
/*
  Purpose:

    HQR computes all eigenvalues of a real upper Hessenberg matrix.

  Discussion:

    HQR finds the eigenvalues of a real
    upper Hessenberg matrix by the QR method.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Martin, Peters, James Wilkinson,
    HQR,
    Numerische Mathematik,
    Volume 14, pages 219-231, 1970.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, two integers determined by 
    BALANC.  If BALANC is not used, set LOW=0, IGH=N-1.

    Input/output, double H[N,N], the N by N upper Hessenberg matrix.
    Information about the transformations used in the reduction to
    Hessenberg form by ELMHES or ORTHES, if performed, is stored
    in the remaining triangle under the Hessenberg matrix.
    On output, the information in H has been destroyed.

    Output, double WR[N], WI[N], the real and imaginary parts of the
    eigenvalues.  The eigenvalues are unordered, except that complex
    conjugate pairs of values appear consecutively, with the eigenvalue
    having positive imaginary part listed first.  If an error break;
    occurred, then the eigenvalues should be correct for indices
    IERR+1 through N.

    Output, int HQR, error flag.
    0, no error.
    J, the limit of 30*N iterations was reached while searching for
      the J-th eigenvalue.
*/
{
  int en;
  int enm2;
  int i;
  int ierr;
  int itn;
  int its;
  int j;
  int k;
  int l;
  int m;
  int na;
  double norm;
  bool notlas;
  double p = 0.0 ;
  double q = 0.0 ;
  double r = 0.0 ;
  double s;
  double t;
  double tst1;
  double tst2;
  double w;
  double x;
  double y;
  double zz;

  ierr = 0;
  norm = 0.0;
  k = 0;
/*
  Store roots isolated by BALANC and compute matrix norm.
*/
  for ( i = 0; i < n; i++ )
  {
    for ( j = k; j < n; j++ )
    {
      norm = norm + fabs ( h[i+j*n] );
    }

    k = i;
    if ( i < low || igh < i )
    {
      wr[i] = h[i+i*n];
      wi[i] = 0.0;
    }
  }

  en = igh;
  t = 0.0;
  itn = 30 * n;
/*
  Search for next eigenvalues.
*/
  if ( igh < low )
  {
    return ierr;
  }

  its = 0;
  na = igh - 1;
  enm2 = igh - 2;
/*
  Look for a single small sub-diagonal element.
*/
  while ( true )
  {
    for ( l = en; low <= l; l-- )
    {
      if ( l == low )
      {
        break;
      }
      s = fabs ( h[l-1+(l-1)*n] ) + fabs ( h[l+l*n] );
      if ( s == 0.0 )
      {
        s = norm;
      }
      tst1 = s;
      tst2 = tst1 + fabs ( h[l+(l-1)*n] );
      if ( tst2 == tst1 )
      {
        break;
      }
    }
/*
  Form shift.
*/
    x = h[en+en*n];
/*
  One root found.
*/
    if ( l == en )
    {
      wr[en] = x + t;
      wi[en] = 0.0;
      en = na;
      if ( en < low )
      {
        return ierr;
      }
      its = 0;
      na = en - 1;
      enm2 = na - 1;
      continue;
    }

    y = h[na+na*n];
    w = h[en+na*n] * h[na+en*n];
/*
  Two roots found.
*/
    if ( l == na )
    {
      p = ( y - x ) / 2.0;
      q = p * p + w;
      zz = sqrt ( fabs ( q ) );
      x = x + t;
/*
  Real root, or complex pair.
*/
      if ( 0.0 <= q )
      {
        zz = p + fabs ( zz ) * r8_sign ( p );
        wr[na] = x + zz;
        if ( zz == 0.0 )
        {
          wr[en] = wr[na];
        }
        else
        {
          wr[en] = x - w / zz;
        }
        wi[na] = 0.0;
        wi[en] = 0.0;
      }
      else
      {
        wr[na] = x + p;
        wr[en] = x + p;
        wi[na] = zz;
        wi[en] = - zz;
      }

      en = enm2;

      if ( en < low )
      {
        return ierr;
      }

      its = 0;
      na = en - 1;
      enm2 = na - 1;
      continue;
    }

    if ( itn == 0 )
    {
      ierr = en;
      return ierr;
    }
/*
  Form an exceptional shift.
*/
    if ( its == 10 || its == 20 )
    {
      t = t + x;

      for ( i = low; i <= en; i++ )
      {
        h[i+i*n] = h[i+i*n] - x;
      }

      s = fabs ( h[en+na*n] ) + fabs ( h[na+enm2*n] );
      x = 0.75 * s;
      y = x;
      w = - 0.4375 * s * s;
    }

    its = its + 1;
    itn = itn - 1;
/*
  Look for two consecutive small sub-diagonal elements.
*/
    for ( m = enm2; l <= m; m-- )
    {
      zz = h[m+m*n];
      r = x - zz;
      s = y - zz;
      p = ( r * s - w ) / h[m+1+m*n] + h[m+(m+1)*n];
      q = h[m+1+(m+1)*n] - zz - r - s;
      r = h[m+2+(m+1)*n];
      s = fabs ( p ) + fabs ( q ) + fabs ( r );
      p = p / s;
      q = q / s;
      r = r / s;

      if ( m == l )
      {
        break;
      }

      tst1 = fabs ( p ) * ( fabs ( h[m-1+(m-1)*n] ) + fabs ( zz ) 
        + fabs ( h[m+1+(m+1)*n] ) );
      tst2 = tst1 + fabs ( h[m+(m-1)*n] ) * ( fabs ( q ) + fabs ( r ) );

      if ( tst2 == tst1 )
      {
        break;
      }
    }

    for ( i = m + 2; i <= en; i++ )
    {
      h[i+(i-2)*n] = 0.0;
      if ( i != m + 2 )
      {
        h[i+(i-3)*n] = 0.0;
      }
    }
/*
  Double QR step involving rows l to EN and columns M to EN.
*/
    for ( k = m; k <= na; k++ )
    {
      notlas = ( k != na );

      if ( k != m )
      {
        p = h[k+(k-1)*n];
        q = h[k+1+(k-1)*n];

        if ( notlas )
        {
          r = h[k+2+(k-1)*n];
        }
        else
        {
          r = 0.0;
        }

        x = fabs ( p ) + fabs ( q ) + fabs ( r );

        if ( x == 0.0 )
        {
          continue;
        }
        p = p / x;
        q = q / x;
        r = r / x;
      }

      s = sqrt ( p * p + q * q + r * r ) * r8_sign ( p );

      if ( k != m )
      {
        h[k+(k-1)*n] = - s * x;
      }
      else if ( l != m )
      {
        h[k+(k-1)*n] = - h[k+(k-1)*n];
      }

      p = p + s;
      x = p / s;
      y = q / s;
      zz = r / s;
      q = q / p;
      r = r / p;
/*
  Row modification.
*/
      if ( ! notlas )
      {
        for ( j = k; j < n; j++ )
        {
          p = h[k+j*n] + q * h[k+1+j*n];
          h[k+j*n] = h[k+j*n] - p * x;
          h[k+1+j*n] = h[k+1+j*n] - p * y;
        }

        j = i4_min ( en, k + 3 );
/*
  Column modification.
*/
        for ( i = 0; i <= j; i++ )
        {
          p = x * h[i+k*n] + y * h[i+(k+1)*n];
          h[i+k*n] = h[i+k*n] - p;
          h[i+(k+1)*n] = h[i+(k+1)*n] - p * q;
        }
      }
/*
  Row modification.
*/
      else
      {
        for ( j = k; j < n; j++ )
        {
          p = h[k+j*n] + q * h[k+1+j*n] + r * h[k+2+j*n];
          h[k+j*n] = h[k+j*n] - p * x;
          h[k+1+j*n] = h[k+1+j*n] - p * y;
          h[k+2+j*n] = h[k+2+j*n] - p * zz;
        }

        j = i4_min ( en, k + 3 );
/*
  Column modification.
*/
        for ( i = 0; i <= j; i++ )
        {
          p = x * h[i+k*n] + y * h[i+(k+1)*n] + zz * h[i+(k+2)*n];
          h[i+k*n] = h[i+k*n] - p;
          h[i+(k+1)*n] = h[i+(k+1)*n] - p * q;
          h[i+(k+2)*n] = h[i+(k+2)*n] - p * r;
        }
      }
    }
  }

  return ierr;
}
/******************************************************************************/

int hqr2 ( int n, int low, int igh, double h[], double wr[], double wi[], 
  double z[] )

/******************************************************************************/
/*
  Purpose:

    HQR2 computes eigenvalues and eigenvectors of a real upper Hessenberg matrix.

  Discussion:

    HQR2 finds the eigenvalues and eigenvectors of a real upper Hessenberg 
    matrix by the QR method.

    The eigenvectors of a real general matrix can also be found
    if ELMHES and ELTRAN or ORTHES and ORTRAN have
    been used to reduce this general matrix to Hessenberg form
    and to accumulate the similarity transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 March 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, determined by the balancing
    routine BALANC.  If BALANC has not been used, set LOW = 0, IGH = N - 1.

    Input/output, double H[N,N], the N by N upper Hessenberg matrix.
    On output, the information in H has been destroyed.

    Output, double WR[N], WI[N], the real and imaginary parts of the
    eigenvalues.  The eigenvalues are unordered, except that complex
    conjugate pairs of values appear consecutively, with the eigenvalue
    having positive imaginary part listed first.  If an error break;
    occurred, then the eigenvalues should be correct for indices
    IERR+1 through N.

    Input/output, double Z[N,N].  On input, the transformation 
    matrix produced by ELTRAN after the reduction by ELMHES, or by ORTRAN after
    the reduction by ORTHES, if performed.  If the eigenvectors of the 
    Hessenberg matrix are desired, Z must contain the identity matrix.  On 
    output, Z contains the real and imaginary parts of the eigenvectors.
    If the I-th eigenvalue is real, the I-th column of Z contains its
    eigenvector.  If the I-th eigenvalue is complex with positive imaginary
    part, the I-th and (I+1)-th columns of Z contain the real and imaginary
    parts of its eigenvector.  The eigenvectors are unnormalized.  If an
    error break; is made, none of the eigenvectors has been found.

    Output, int HQR2, error flag.
    0, for normal return,
    J, if the limit of 30*N iterations is exhausted while the J-th
      eigenvalue is being sought.
*/
{
  int en;
  int enm2;
  int i;
  int ierr;
  int itn;
  int its;
  int j;
  int k;
  int l;
  int m;
  int na;
  double norm;
  bool notlas;
  double p = 0.0 ;
  double q = 0.0;
  double r = 0.0;
  double ra;
  double s = 0.0;
  double sa;
  double t;
  double ti;
  double tr;
  double tst1;
  double tst2;
  double vi;
  double vr;
  double w;
  double x;
  double y;
  double zz = 0.0 ;

  ierr = 0;
  norm = 0.0;
  k = 0;
/*
  Store roots isolated by BALANC and compute the matrix norm.
*/
  for ( i = 0; i < n; i++ )
  {
    for ( j = k; j < n; j++ )
    {
      norm = norm + fabs ( h[i+j*n] );
    }
    k = i;
    if ( i < low || igh < i )
    {
      wr[i] = h[i+i*n];
      wi[i] = 0.0;
    }
  }

  en = igh;
  t = 0.0;
  itn = 30 * n;
/*
  Search for next eigenvalues.
*/
  while ( low <= en )
  {
    its = 0;
    na = en - 1;
    enm2 = na - 1;
/*
  Look for single small sub-diagonal element.
*/
   while ( true )
   {
      for ( l = en; low <= l; l-- )
      {
        if ( l == low )
        {
          break;
        }
        s = fabs ( h[l-1+(l-1)*n] ) + fabs ( h[l+l*n] );
        if ( s == 0.0 )
        {
          s = norm;
        }
        tst1 = s;
        tst2 = tst1 + fabs ( h[l+(l-1)*n] );
        if ( tst2 == tst1 )
        {
          break;
        }
      }
/*
  Form shift.
*/
      x = h[en+en*n];
/*
  One root found.
*/
      if ( l == en )
      {
        h[en+en*n] = x + t;
        wr[en] = h[en+en*n];
        wi[en] = 0.0;
        en = na;
        break;
      }

      y = h[na+na*n];
      w = h[en+na*n] * h[na+en*n];
/*
  Two roots found.
*/
      if ( l == na )
      {
        p = ( y - x ) / 2.0;
        q = p * p + w;
        zz = sqrt ( fabs ( q ) );
        h[en+en*n] = x + t;
        x = h[en+en*n];
        h[na+na*n] = y + t;
/*
  Complex pair.
*/
        if ( q < 0.0 )
        {
          wr[na] = x + p;
          wr[en] = x + p;
          wi[na] = zz;
          wi[en] = - zz;
        }
/*
  Real pair.
*/
        else
        {
          zz = p + fabs ( zz ) * r8_sign ( p );
          wr[na] = x + zz;
          wr[en] = wr[na];

          if ( zz != 0.0 )
          {
            wr[en] = x - w / zz;
          }

          wi[na] = 0.0;
          wi[en] = 0.0;
          x = h[en+na*n];
          s = fabs ( x ) + fabs ( zz );
          p = x / s;
          q = zz / s;
          r = sqrt ( p * p + q * q );
          p = p / r;
          q = q / r;
/*
  Row modification.
*/
          for ( j = na; j < n; j++ )
          {
            zz = h[na+j*n];
            h[na+j*n] = q * zz + p * h[en+j*n];
            h[en+j*n] = q * h[en+j*n] - p * zz;
          }
/*
  Column modification.
*/
          for ( i = 0; i <= en; i++ )
          {
            zz = h[i+na*n];
            h[i+na*n] = q * zz + p * h[i+en*n];
            h[i+en*n] = q * h[i+en*n] - p * zz;
          }
/*
  Accumulate transformations.
*/
         for ( i = low; i <= igh; i++ )
         {
            zz = z[i+na*n];
            z[i+na*n] = q * zz + p * z[i+en*n];
            z[i+en*n] = q * z[i+en*n] - p * zz;
          }
        }
        en = enm2;
        break;
      }

      if ( itn == 0 )
      {
        ierr = en;
        return ierr;
      }
/*
  Form exceptional shift.
*/
      if ( its == 10 || its == 20 )
      {
        t = t + x;
        for ( i = low; i <= en; i++ )
        {
          h[i+i*n] = h[i+i*n] - x;
        }
        s = fabs ( h[en+na*n] ) + fabs ( h[na+enm2*n] );
        x = 0.75 * s;
        y = x;
        w = - 0.4375 * s * s;
      }

      its = its + 1;
      itn = itn - 1;
/*
  Look for two consecutive small sub-diagonal elements.
*/
      for ( m = enm2; l <= m; m-- )
      {
        zz = h[m+m*n];
        r = x - zz;
        s = y - zz;
        p = ( r * s - w ) / h[m+1+m*n] + h[m+(m+1)*n];
        q = h[m+1+(m+1)*n] - zz - r - s;
        r = h[m+2+(m+1)*n];
        s = fabs ( p ) + fabs ( q ) + fabs ( r );
        p = p / s;
        q = q / s;
        r = r / s;
        if ( m == l )
        {
          break;
        }
        tst1 = fabs ( p ) * ( fabs ( h[m-1+(m-1)*n] ) + fabs ( zz ) 
          + fabs ( h[m+1+(m+1)*n] ) );
        tst2 = tst1 + fabs ( h[m+(m-1)*n] ) * ( fabs ( q ) + fabs ( r ) );
        if ( tst2 == tst1 )
        {
          break;
        }
      }

      for ( i = m + 2; i <= en; i++ )
      {
        h[i+(i-2)*n] = 0.0;
        if ( i != m + 2 )
        {
          h[i+(i-3)*n] = 0.0;
        }
      }
/*
  Double QR step involving rows L to EN and columns M to EN.
*/
      for ( k = m; k <= na; k++ )
      {
        notlas = ( k != na );

        if ( k != m )
        {
          p = h[k+(k-1)*n];
          q = h[k+1+(k-1)*n];
          r = 0.0;
          if ( notlas )
          {
            r = h[k+2+(k-1)*n];
          }
          x = fabs ( p ) + fabs ( q ) + fabs ( r );
          if ( x == 0.0 )
          {
            continue;
          }
          p = p / x;
          q = q / x;
          r = r / x;
        }

        s = sqrt ( p * p + q * q + r * r ) * r8_sign ( p );

        if ( k != m )
        {
          h[k+(k-1)*n] = - s * x;
        }
        else if ( l != m )
        {
          h[k+(k-1)*n] = - h[k+(k-1)*n];
        }

        p = p + s;
        x = p / s;
        y = q / s;
        zz = r / s;
        q = q / p;
        r = r / p;

        if ( ! notlas )
        {
/*
  Row modification.
*/
          for ( j = k; j < n; j++ )
          {
            p = h[k+j*n] + q * h[k+1+j*n];
            h[k+j*n] = h[k+j*n] - p * x;
            h[k+1+j*n] = h[k+1+j*n] - p * y;
          }

          j = i4_min ( en, k + 3 );
/*
  Column modification.
*/
          for ( i = 0; i <= j; i++ )
          {
            p = x * h[i+k*n] + y * h[i+(k+1)*n];
            h[i+k*n] = h[i+k*n] - p;
            h[i+(k+1)*n] = h[i+(k+1)*n] - p * q;
          }
/*
  Accumulate transformations.
*/
          for ( i = low; i <= igh; i++ )
          {
            p = x * z[i+k*n] + y * z[i+(k+1)*n];
            z[i+k*n] = z[i+k*n] - p;
            z[i+(k+1)*n] = z[i+(k+1)*n] - p * q;
          }
        }
/*
  Row modification.
*/
        else
        {
          for ( j = k; j < n; j++ )
          {
            p = h[k+j*n] + q * h[k+1+j*n] + r * h[k+2+j*n];
            h[k+j*n] = h[k+j*n] - p * x;
            h[k+1+j*n] = h[k+1+j*n] - p * y;
            h[k+2+j*n] = h[k+2+j*n] - p * zz;
          }

          j = i4_min ( en, k + 3 );
/*
  Column modification.
*/
          for ( i = 0; i <= j; i++ )
          {
            p = x * h[i+k*n] + y * h[i+(k+1)*n] + zz * h[i+(k+2)*n];
            h[i+k*n] = h[i+k*n] - p;
            h[i+(k+1)*n] = h[i+(k+1)*n] - p * q;
            h[i+(k+2)*n] = h[i+(k+2)*n] - p * r;
          }
/*
  Accumulate transformations.
*/
          for ( i = low; i <= igh; i++ )
          {
            p = x * z[i+k*n] + y * z[i+(k+1)*n] + zz * z[i+(k+2)*n];
            z[i+k*n] = z[i+k*n] - p;
            z[i+(k+1)*n] = z[i+(k+1)*n] - p * q;
            z[i+(k+2)*n] = z[i+(k+2)*n] - p * r;
          }
        }
      }
    }
  }
/*
  All roots found.
  Backsubstitute to find vectors of upper triangular form.
*/
  if ( norm == 0.0 )
  {
    return ierr;
  }

  for ( en = n - 1; 0 <= en; en-- )
  {
    p = wr[en];
    q = wi[en];
    na = en - 1;

    if ( 0.0 < q )
    {
      continue;
    }
/*
  Real vector
*/
    else if ( q == 0.0 )
    {
      m = en;
      h[en+en*n] = 1.0;

      for ( i = en - 1; en - na - 1 <= i; i-- )
      {
        w = h[i+i*n] - p;
        r = 0.0;
        for ( j = m; j <= en; j++ )
        {
          r = r + h[i+j*n] * h[j+en*n];
        }

        if ( wi[i] < 0.0 )
        {
          zz = w;
          s = r;
          continue;
        }

        m = i;

        if ( wi[i] == 0.0 )
        {
          t = w;

          if ( t == 0.0 )
          {
            tst1 = norm;
            t = tst1;

            while ( true )
            {
              t = 0.01 * t;
              tst2 = norm + t;
              if ( tst2 <= tst1 )
              {
                break;
              }
            }
          }
          h[i+en*n] = - r / t;
        }
/*
  Solve real equations.
*/
        else
        {
          x = h[i+(i+1)*n];
          y = h[i+1+i*n];
          q = ( wr[i] - p ) * ( wr[i] - p ) + wi[i] * wi[i];
          t = ( x * s - zz * r ) / q;
          h[i+en*n] = t;

          if ( fabs ( zz ) < fabs ( x ) )
          {
            h[i+1+en*n] = ( - r - w * t ) / x;
          }
          else
          {
            h[i+1+en*n] = ( - s - y * t ) / zz;
          }
        }
/*
  Overflow control.
*/
        t = fabs ( h[i+en*n] );
        if ( t != 0.0 )
        {
          tst1 = t;
          tst2 = tst1 + 1.0 / tst1;
          if ( tst2 <= tst1 )
          {
            for ( j = i; j <= en; j++ )
            {
              h[j+en*n] = h[j+en*n] / t;
            }
          }
        }
      }
    }
/*
  Complex vector
*/
    else if ( q < 0.0 )
    {
      m = na;
/*
  Last vector component chosen imaginary, so that the eigenvector
  matrix is triangular.
*/
      if ( fabs ( h[na+en*n] ) < fabs ( h[en+na*n] ) )
      {
        h[na+na*n] = q / h[en+na*n];
        h[na+en*n] = - ( h[en+en*n] - p ) / h[en+na*n];
      }
      else
      {
        cdiv ( 0.0, -h[na+en*n], h[na+na*n] - p, q, &tr, &ti );
        h[na+na*n] = tr;
        h[na+en*n] = ti;
      }

      h[en+na*n] = 0.0;
      h[en+en*n] = 1.0;
      enm2 = na - 1;

      for ( i = na - 1; na - enm2 <= i; i-- )
      {
        w = h[i+i*n] - p;
        ra = 0.0;
        sa = 0.0;
        for ( j = m; j <= en; j++ )
        {
          ra = ra + h[i+j*n] * h[j+na*n];
          sa = sa + h[i+j*n] * h[j+en*n];
        }
        if ( wi[i] < 0.0 )
        {
          zz = w;
          r = ra;
          s = sa;
        }

        m = i;

        if ( wi[i] == 0.0 )
        {
          cdiv ( -ra, -sa, w, q, &tr, &ti );
          h[i+na*n] = tr;
          h[i+en*n] = ti;
        }
/*
  Solve complex equations.
*/
        else
        {
          x = h[i+(i+1)*n];
          y = h[i+1+i*n];
          vr = ( wr[i] - p ) * ( wr[i] - p ) + wi[i] * wi[i] - q * q;
          vi = ( wr[i] - p ) * 2.0 * q;

          if ( vr == 0.0 && vi == 0.0 )
          {
            tst1 = norm * ( fabs ( w ) + fabs ( q ) + fabs ( x ) 
              + fabs ( y ) + fabs ( zz ) );
            vr = tst1;

            while ( true )
            {
              vr = 0.01 * vr;
              tst2 = tst1 + vr;
              if ( tst2 <= tst1 )
              {
                break;
              }
            }
          }

          cdiv ( x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, 
            vr, vi, &tr, &tr );
          h[i+na*n] = tr;
          h[i+en*n] = ti;

          if ( fabs ( zz ) + fabs ( q ) < fabs ( x ) )
          {
            h[i+1+na*n] = ( - ra - w * h[i+na*n] + q * h[i+en*n] ) / x;
            h[i+1+en*n] = ( - sa - w * h[i+en*n] - q * h[i+na*n] ) / x;
          }
          else
          {
            cdiv ( - r - y * h[i+na*n], - s - y * h[i+en*n], zz, q, 
              &tr, &ti );
            h[i+1+na*n] = tr;
            h[i+1+en*n] = ti;
          }
        }
/*
  Overflow control.
*/
        t = r8_max ( fabs ( h[i+na*n] ), fabs ( h[i+en*n] ) );

        if ( t != 0.0 )
        {
          tst1 = t;
          tst2 = tst1 + 1.0 / tst1;
          if ( tst2 <= tst1 )
          {
            for ( j = i; j <= en; j++ )
            {
              h[j+na*n] = h[j+na*n] / t;
              h[j+en*n] = h[j+en*n] / t;
            }
          }
        }
      }
    }
  }
/*
  End back substitution.

  Vectors of isolated roots.
*/
  for ( i = 0; i < n; i++ )
  {
    if ( i < low || igh < i )
    {
      for ( j = i; j < n; j++ )
      {
        z[i+j*n] = h[i+j*n];
      }
    }
  }
/*
  Multiply by transformation matrix to give vectors of original full matrix.
*/
  for ( j = n - 1; low <= j; j-- )
  {
    m = i4_min ( j, igh );
    for ( i = low; i <= igh; i++ )
    {
      zz = 0.0;
      for ( k = low; k <= m; k++ )
      {
        zz = zz + z[i+k*n] * h[k+j*n];
      }
      z[i+j*n] = zz;
    }
  }

  return ierr;
}
/******************************************************************************/

void htrib3 ( int n, double a[], double tau[], int m, double zr[], double zi[] )

/******************************************************************************/
/*
  Purpose:

    HTRIB3 determines eigenvectors by undoing the HTRID3 transformation.

  Discussion:

    HTRIB3 forms the eigenvectors of a complex hermitian
    matrix by back transforming those of the corresponding
    real symmetric tridiagonal matrix determined by HTRID3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, is the order of the matrix.

    Input, double A[N,N], contains information about the unitary
    transformations used in the reduction by HTRID3.

    Input, double TAU[2,N], contains further information about the
    transformations.

    Input, int M, the number of eigenvectors to be back
    transformed.

    Input/output, double ZR[N,M], ZI[N,M].  On input, ZR contains 
    the eigenvectors to be back transformed.  On output, ZR and ZI contain
    the real and imaginary parts of the transformed eigenvectors.
*/
{
  double h;
  int i;
  int j;
  int k;
  double s;
  double si;

  if ( m == 0 )
  {
    return;
  }
/*
  Transform the eigenvectors of the real symmetric tridiagonal matrix
  to those of the hermitian tridiagonal matrix.
*/
  for ( k = 0; k < n; k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      zi[k+j*n] = - zr[k+j*n] * tau[1+k*2];
      zr[k+j*n] =   zr[k+j*n] * tau[0+k*2];
    }
  }
/*
  Recover and apply the Householder matrices.
*/
  for ( i = 1; i < n; i++ )
  {
    h = a[i+i*n];

    if ( h != 0.0 )
    {
      for ( j = 0; j < m; j++ )
      {
        s = 0.0;
        si = 0.0;

        for ( k = 0; k < i; k++ )
        {
          s  = s  + a[i+k*n] * zr[k+j*n] - a[k+i*n] * zi[k+j*n];
          si = si + a[i+k*n] * zi[k+j*n] + a[k+i*n] * zr[k+j*n];
        }

        s = ( s / h ) / h;
        si = ( si / h ) / h;

        for ( k = 0; k < i; k++ )
        {
          zr[k+j*n] = zr[k+j*n] - s  * a[i+k*n] - si * a[k+i*n];
          zi[k+j*n] = zi[k+j*n] - si * a[i+k*n] + s  * a[k+i*n];
        }
      }
    }
  }

  return;
}
/******************************************************************************/

void htribk ( int n, double ar[], double ai[], double tau[], int m, double zr[], 
  double zi[] )

/******************************************************************************/
/*
  Purpose:

    HTRIBK determines eigenvectors by undoing the HTRIDI transformation.

  Discussion:

    HTRIBK forms the eigenvectors of a complex hermitian
    matrix by back transforming those of the corresponding
    real symmetric tridiagonal matrix determined by HTRIDI.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double AR[N,N], AI[N,N], contain information about
    the unitary transformations used in the reduction by HTRIDI in their
    full lower triangles, except for the diagonal of AR.

    Input, double TAU[2,N], contains further information about the
    transformations.

    Input, int  M, the number of eigenvectors to be back
    transformed.

    Input/output, double ZR[N,M], ZI[N,M].  On input, ZR contains 
    the eigenvectors to be back transformed.  On output, ZR and ZI contain
    the real and imaginary parts of the transformed eigenvectors.
*/
{
  double h;
  int i;
  int j;
  int k;
  double s;
  double si;

  if ( m == 0 )
  {
    return;
  }
/*
  Transform the eigenvectors of the real symmetric tridiagonal matrix to
  those of the hermitian tridiagonal matrix.
*/
  for ( k = 0; k < n; k++ )
  {
    for ( j = 0; j < m; j++ )
    {
      zi[k+j*n] = - zr[k+j*n] * tau[1+k*2];
      zr[k+j*n] =   zr[k+j*n] * tau[0+k*2];
    }
  }
/*
  Recover and apply the Householder matrices.
*/
  for ( i = 1; i < n; i++ )
  {
    h = ai[i+i*n];

    if ( h != 0.0 )
    {
      for ( j = 0; j < m; j++ )
      {
        s = 0.0;
        si = 0.0;
        for ( k = 0; k < i; k++ )
        {
          s =  s  + ar[i+k*n] * zr[k+j*n] - ai[i+k*n] * zi[k+j*n];
          si = si + ar[i+k*n] * zi[k+j*n] + ai[i+k*n] * zr[k+j*n];
        }

        s = ( s / h ) / h;
        si = ( si / h ) / h;

        for ( k = 0; k < i; k++ )
        {
          zr[k+j*n] = zr[k+j*n] - s  * ar[i+k*n] - si * ai[i+k*n];
          zi[k+j*n] = zi[k+j*n] - si * ar[i+k*n] + s  * ai[i+k*n];
        }
      }
    }
  }

  return;
}
/******************************************************************************/

void htrid3 ( int n, double a[], double d[], double e[], double e2[], 
  double tau[] )

/******************************************************************************/
/*
  Purpose:

    HTRID3 tridiagonalizes a complex hermitian packed matrix.

  Discussion:

    HTRID3 reduces a complex hermitian matrix, stored as a single square 
    array, to a real symmetric tridiagonal matrix using unitary similarity 
    transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double A[N,N].  On input, the lower triangle of 
    the complex hermitian input matrix.  The real parts of the matrix elements 
    are stored in the full lower triangle of A, and the imaginary parts are 
    stored in the transposed positions of the strict upper triangle of A.  No 
    storage is required for the zero imaginary parts of the diagonal elements.
    On output, A contains information about the unitary transformations
    used in the reduction.

    Output, double D[N], the diagonal elements of the
    tridiagonal matrix.

    Output, double E[N], the subdiagonal elements of the tridiagonal
    matrix in E(2:N).  E(1) is set to zero.

    Output, double E2[N], the squares of the corresponding elements 
    of E.  E2 may coincide with E if the squares are not needed.

    Output, double TAU[2,N], contains further information about the
    transformations.
*/
{
  double f;
  double fi;
  double g;
  double gi;
  double h;
  double hh;
  int i;
  int j;
  int k;
  double scale;
  double si;

  tau[0+(n-1)*2] = 1.0;
  tau[1+(n-1)*2] = 0.0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    h = 0.0;
    scale = 0.0;

    if ( i < 1 )
    {
      e[i] = 0.0;
      e2[i] = 0.0;
    }
    else
    {
/*
  Scale row.
*/
      for ( k = 0; k < i; k++ )
      {
        scale = scale + fabs ( a[i+k*n] ) + fabs ( a[k+i*n] );
      }

      if ( scale == 0.0 )
      {
        tau[0+(i-1)*2] = 1.0;
        tau[1+(i-1)*2] = 0.0;
        e[i] = 0.0;
        e2[i] = 0.0;
      }
      else
      {
        for ( k = 0; k < i; k++ )
        {
          a[i+k*n] = a[i+k*n] / scale;
          a[k+i*n] = a[k+i*n] / scale;
          h = h + a[i+k*n] * a[i+k*n] + a[k+i*n] * a[k+i*n];
        }

        e2[i] = scale * scale * h;
        g = sqrt ( h );
        e[i] = scale * g;
        f = pythag ( a[i+(i-1)*n], a[i-1+i*n] );
/*
  Form next diagonal element of matrix T.
*/
        if ( f != 0.0 )
        {
          tau[0+(i-1)*2] = ( a[i-1+i*n]   * tau[1+i*2] - a[i+(i-1)*n] * tau[0+i*2] ) / f;
          si =             ( a[i+(i-1)*n] * tau[1+i*2] + a[i-1+i*n]   * tau[0+i*2] ) / f;
          h = h + f * g;
          g = 1.0 + g / f;
          a[i+(i-1)*n] = g * a[i+(i-1)*n];
          a[i-1+i*n]   = g * a[i-1+i*n];

          if ( i == 1 )
          {
            a[1+0*n] = scale * a[1+0*n];
            a[0+1*n] = scale * a[0+1*n];
            tau[1+(i-1)*2] = - si;
            d[i] = a[i+i*n];
            a[i+i*n] = scale * sqrt ( h );
            continue;
          }
        }
        else
        {
          tau[0+(i-1)*2] = - tau[0+i*2];
          si = tau[1+i*2];
          a[i+(i-1)*n] = g;
        }

        f = 0.0;

        for ( j = 0; j < i; j++ )
        {
          g = 0.0;
          gi = 0.0;
/*
  Form element of A*U.
*/
          for ( k = 0; k < j; k++ )
          {
            g  = g  + a[j+k*n] * a[i+k*n] + a[k+j*n] * a[k+i*n];
            gi = gi - a[j+k*n] * a[k+i*n] + a[k+j*n] * a[i+k*n];
          }

          g  = g  + a[j+j*n] * a[i+j*n];
          gi = gi - a[j+j*n] * a[j+i*n];

          for ( k = j + 1; k < i; k++ )
          {
            g  = g  + a[k+j*n] * a[i+k*n] - a[j+k*n] * a[k+i*n];
            gi = gi - a[k+j*n] * a[k+i*n] - a[j+k*n] * a[i+k*n];
          }
/*
  Form element of P.
*/
          e[j] = g / h;
          tau[1+j*2] = gi / h;
          f = f + e[j] * a[i+j*n] - tau[1+j*2] * a[j+i*n];
        }

        hh = f / ( h + h );
/*
  Form reduced A.
*/
        for ( j = 0; j < i; j++ )
        {
          f = a[i+j*n];
          g = e[j] - hh * f;
          e[j] = g;
          fi = - a[j+i*n];
          gi = tau[1+j*2] - hh * fi;
          tau[1+j*2] = - gi;
          a[j+j*n] = a[j+j*n] - 2.0 * ( f * g + fi * gi );

          for ( k = 0; k < j; k++ )
          {
            a[j+k*n] = a[j+k*n] 
              - f * e[k] - g * a[i+k*n] + fi * tau[1+k*2] + gi * a[k+i*n];
            a[k+j*n] = a[k+j*n] 
              - f * tau[1+k*2] - g * a[k+i*n] - fi * e[k] - gi * a[i+k*n];
          }
        }

        for ( j = 0; j < i; j++ )
        {
          a[i+j*n] = scale * a[i+j*n];
          a[j+i*n] = scale * a[j+i*n];
        }
        tau[1+(i-1)*2] = - si;
      }
    }

    d[i] = a[i+i*n];
    a[i+i*n] = scale * sqrt ( h );
  }

  return;
}
/******************************************************************************/

void htridi ( int n, double ar[], double ai[], double d[], double e[], 
  double e2[], double tau[] )

/******************************************************************************/
/*
  Purpose:

    HTRIDI tridiagonalizes a complex hermitian matrix.

  Discussion:

    HTRIDI reduces a complex hermitian matrix to a real symmetric
    tridiagonal matrix using unitary similarity transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double AR[N,N], AI[N,N].  On input, the real
    and imaginary parts, respectively, of the complex hermitian input matrix.
    Only the lower triangle of the matrix need be supplied.
    On output, information about the unitary transformations used in the
    reduction in their full lower triangles.  Their strict upper triangles
    and the diagonal of AR are unaltered.

    Output, double D[N], the diagonal elements of the
    tridiagonal matrix.

    Output, double E[N], the subdiagonal elements of the tridiagonal
    matrix in its last N-1 positions.  E(1) is set to zero.

    Output, double E2[N], the squares of the corresponding elements 
    of E.  E2 may coincide with E if the squares are not needed.

    Output, double TAU[2,N], contains further information about the
    transformations.
*/
{
  double f;
  double fi;
  double g;
  double gi;
  double h;
  double hh;
  int i;
  int j;
  int k;
  int l;
  double scale;
  double si;

  tau[0+(n-1)*2] = 1.0;
  tau[1+(n-1)*2] = 0.0;

  for ( i = 0; i < n; i++ )
  {
    d[i] = ar[i+i*n];
  }

  for ( i = n - 1; 0 <= i; i-- )
  {
    l = i - 1;
    h = 0.0;
    scale = 0.0;

    if ( i == 0 )
    {
      e[i] = 0.0;
      e2[i] = 0.0;
      hh = d[i];
      d[i] = ar[i+i*n];
      ar[i+i*n] = hh;
      ai[i+i*n] = scale * sqrt ( h );
    }
/*
  Scale row.
*/
    else
    {
      for ( k = 0; k < i; k++ )
      {
        scale = scale + fabs ( ar[i+k*n] ) + fabs ( ai[i+k*n] );
      }

      if ( scale == 0.0 )
      {
        tau[0+(i-1)*2] = 1.0;
        tau[1+(i-1)*2] = 0.0;
        e[i] = 0.0;
        e2[i] = 0.0;
        hh = d[i];
        d[i] = ar[i+i*n];
        ar[i+i*n] = hh;
        ai[i+i*n] = scale * sqrt ( h );
        continue;
      }

      for ( j = 0; j < i; j++ )
      {
        ar[i+j*n] = ar[i+j*n] / scale;
        ai[i+j*n] = ai[i+j*n] / scale;
      }

      for ( k = 0; k < i; k++ )
      {
        h = h + ar[i+k*n] * ar[i+k*n] + ai[i+k*n] * ai[i+k*n];
      }

      e2[i] = scale * scale * h;
      g = sqrt ( h );
      e[i] = scale * g;
      f = pythag ( ar[i+(i-1)*n], ai[i+(i-1)*n] );
/*
  Form next diagonal element of matrix T.
*/
      if ( f != 0.0 )
      {
        tau[0+(i-1)*2] = ( ai[i+(i-1)*n] * tau[1+i*2] 
                         - ar[i+(i-1)*n] * tau[0+i*2] ) / f;
        si =             ( ar[i+(i-1)*n] * tau[1+i*2] 
                         + ai[i+(i-1)*n] * tau[0+i*2] ) / f;
        h = h + f * g;
        g = 1.0 + g / f;
        ar[i+(i-1)*n] = g * ar[i+(i-1)*n];
        ai[i+(i-1)*n] = g * ai[i+(i-1)*n];

        if ( i == 1 )
        {
          for ( j = 0; j < i; j++ )
          {
            ar[i+j*n] = scale * ar[i+j*n];
            ai[i+j*n] = scale * ai[i+j*n];
          }
          tau[1+(i-1)*2] = - si;
          hh = d[i];
          d[i] = ar[i+i*n];
          ar[i+i*n] = hh;
          ai[i+i*n] = scale * sqrt ( h );
          continue;
        }
      }
      else
      {
        tau[0+(i-1)*2] = - tau[0+i*2];
        si = tau[1+i*2];
        ar[i+(i-1)*n] = g;
      }

      f = 0.0;

      for ( j = 0; j < i; j++ )
      {
        g = 0.0;
        gi = 0.0;
/*
  Form element of A*U.
*/
        for ( k = 0; k <= j; k++ )
        {
          g  = g  + ar[j+k*n] * ar[i+k*n] + ai[j+k*n] * ai[i+k*n];
          gi = gi - ar[j+k*n] * ai[i+k*n] + ai[j+k*n] * ar[i+k*n];
        }

        for ( k = j + 1; k < i; k++ )
        {
          g  = g  + ar[k+j*n] * ar[i+k*n] - ai[k+j*n] * ai[i+k*n];
          gi = gi - ar[k+j*n] * ai[i+k*n] - ai[k+j*n] * ar[i+k*n];
        }
/*
  Form element of P.
*/
        e[j] = g / h;
        tau[1+j*2] = gi / h;
        f = f + e[j] * ar[i+j*n] - tau[1+j*2] * ai[i+j*n];
      }

      hh = f / ( h + h );
/*
  Form the reduced A.
*/
      for ( j = 0; j < i; j++ )
      {
        f = ar[i+j*n];
        g = e[j] - hh * f;
        e[j] = g;
        fi = - ai[i+j*n];
        gi = tau[1+j*2] - hh * fi;
        tau[1+j*2] = - gi;

        for ( k = 0; k <= j; k++ )
        {
          ar[j+k*n] = ar[j+k*n] - f * e[k] - g * ar[i+k*n] + fi * tau[1+k*2]
            + gi * ai[i+k*n];
          ai[j+k*n] = ai[j+k*n] - f * tau[1+k*2] - g * ai[i+k*n] - fi * e[k] 
            - gi * ar[i+k*n];
        }
      }

      for ( j = 0; j < i; j++ )
      {
        ar[i+j*n] = scale * ar[i+j*n];
        ai[i+j*n] = scale * ai[i+j*n];
      }
      tau[1+l*2] = - si;

      hh = d[i];
      d[i] = ar[i+i*n];
      ar[i+i*n] = hh;
      ai[i+i*n] = scale * sqrt ( h );
    }
  }

  return;
}
/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

void i4vec_print ( int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PRINT prints an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 November 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, int A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %6d: %8d\n", i, a[i] );
  }
  return;
}
/******************************************************************************/

int imtql1 ( int n, double d[], double e[] )

/******************************************************************************/
/*
  Purpose:

    IMTQL1 computes all eigenvalues of a symmetric tridiagonal matrix.

  Discussion:

    This routine finds the eigenvalues of a symmetric
    tridiagonal matrix by the implicit QL method.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D[N].  On input, the diagonal elements of
    the matrix.  On output, the eigenvalues in ascending order.  If an error
    exit is made, the eigenvalues are correct and ordered for indices
    1,2,...IERR-1, but may not be the smallest eigenvalues.

    Input/output, double E[N].  On input, the subdiagonal elements
    of the matrix in its last N-1 positions.  E(1) is arbitrary.  On output,
    E has been overwritten.

    Output, int IMTQL1, error flag.
    0, normal return,
    J, if the J-th eigenvalue has not been determined after 30 iterations.
*/
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ierr;
  int its;
  int l;
  int m;
  double p;
  double r;
  double s;
  bool skip;
  double tst1;
  double tst2;

  ierr = 0;

  if ( n == 1 )
  {
    return ierr;
  }

  for ( i = 1; i < n; i++ )
  {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;

  for ( l = 0; l < n; l++ )
  {
    its = 0;
/*
  Look for a small sub-diagonal element.
*/
    while ( true )
    {
      m = l;

      for ( m = l; m < n - 1; m++ )
      {
        tst1 = fabs ( d[m] ) + fabs ( d[m+1] );
        tst2 = tst1 + fabs ( e[m] );

        if ( tst2 == tst1 )
        {
          break;
        }
      }
/*
  Order the eigenvalues.
*/
      p = d[l];

      if ( m == l )
      {
        for ( i = l; 0 <= i; i-- )
        {
          if ( i == 0 )
          {
            d[i] = p;
            break;
          }

          if ( d[i-1] <= p )
          {
            d[i] = p;
            break;
          }
         d[i] = d[i-1];
        }
        break;
      }
      else
      {
        if ( 30 <= its )
        {
          ierr = l + 1;
          return ierr;
        }
        its = its + 1;
/*
  Form shift.
*/
        g = ( d[l+1] - p ) / ( 2.0 * e[l] );
        r = pythag ( g, 1.0 );
        g = d[m] - p + e[l] / ( g + fabs ( r ) * r8_sign ( g ) );
        s = 1.0;
        c = 1.0;
        p = 0.0;

        skip = false;

        for ( i = m - 1; l <= i; i-- )
        {
          f = s * e[i];
          b = c * e[i];
          r = pythag ( f, g );
          e[i+1] = r;
/*
  Recover from underflow.
*/
          if ( r == 0.0 )
          {
            d[i+1] = d[i+1] - p;
            e[m] = 0.0;
            skip = true;
            break;
          }

          s = f / r;
          c = g / r;
          g = d[i+1] - p;
          r = ( d[i] - g ) * s + 2.0 * c * b;
          p = s * r;
          d[i+1] = g + p;
          g = c * r - b;
        }

        if ( ! skip )
        {
          d[l] = d[l] - p;
          e[l] = g;
          e[m] = 0.0;
        }
      }
    }
  }

  return ierr;
}
/******************************************************************************/

int imtql2 ( int n, double d[], double e[], double z[] )

/******************************************************************************/
/*
  Purpose:

    IMTQL2 computes all eigenvalues/vectors of a symmetric tridiagonal matrix.

  Discussion:

    IMTQL2 finds the eigenvalues and eigenvectors of a symmetric tridiagonal 
    matrix by the implicit QL method.

    The eigenvectors of a full symmetric matrix can also be found if TRED2 
    has been used to reduce this full matrix to tridiagonal form.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D[N].  On input, the diagonal elements of
    the input matrix.  On output, the eigenvalues in ascending order.  If an
    error exit is made, the eigenvalues are correct but
    unordered for indices 1,2,...,IERR-1.

    Input/output, double E[N].  On input, the subdiagonal elements
    of the input matrix in E(2:N).  E(1) is arbitrary.  On output, E is
    overwritten.

    Input/output, double Z[N,N].  On input, the transformation
    matrix produced in the reduction by TRED2, if performed.  If the
    eigenvectors of the tridiagonal matrix are desired, Z must contain the
    identity matrix.  On output, Z contains orthonormal eigenvectors of the
    symmetric tridiagonal (or full) matrix.  If an error exit is made, Z
    contains the eigenvectors associated with the stored eigenvalues.

    Output, int IMTQL2, error flag.
    0, for normal return,
    J, if the J-th eigenvalue has not been determined after 30 iterations.
*/
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ierr;
  int ii;
  int its;
  int j;
  int k;
  int l;
  int m;
  double p;
  double r;
  double s;
  double t;
  double tst1;
  double tst2;

  ierr = 0;

  if ( n == 1 )
  {
    return ierr;
  }

  for ( i = 1; i < n; i++ )
  {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;

  for ( l = 0; l < n; l++ )
  {
    its = 0;
/*
  Look for a small sub-diagonal element.
*/
    while ( true )
    {
      m = l;

      for ( m = l; m < n - 1; m++ )
      {
        tst1 = fabs ( d[m] ) + fabs ( d[m+1] );
        tst2 = tst1 + fabs ( e[m] );

        if ( tst2 == tst1 )
        {
          break;
        }
      }

      p = d[l];

      if ( m == l )
      {
        break;
      }

      if ( 30 <= its )
      {
        ierr = l + 1;
        return ierr;
      }

      its = its + 1;
/*
  Form shift.
*/
      g = ( d[l+1] - p ) / ( 2.0 * e[l] );
      r = pythag ( g, 1.0 );
      g = d[m] - p + e[l] / ( g + fabs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;

      for ( i = m - 1; l <= i; i-- )
      {
        f = s * e[i];
        b = c * e[i];
        r = pythag ( f, g );
        e[i+1] = r;
/*
  Recover from underflow.
*/
        if ( r == 0.0 )
        {
          d[i+1] = d[i+1] - p;
          e[m] = 0.0;
          continue;
        }

        s = f / r;
        c = g / r;
        g = d[i+1] - p;
        r = ( d[i] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i+1] = g + p;
        g = c * r - b;
/*
  Form vector.
*/
        for ( k = 0; k < n; k++ )
        {
          f = z[k+(i+1)*n];
          z[k+(i+1)*n] = s * z[k+i*n] + c * f;
          z[k+i*n]     = c * z[k+i*n] - s * f;
        }
      }

      d[l] = d[l] - p;
      e[l] = g;
      e[m] = 0.0;
    }
  }
/*
  Order eigenvalues and eigenvectors.
*/
  for ( i = 0; i < n - 1; i++ )
  {
    k = i;
    p = d[i];

    for ( j = i + 1; j < n; j++ )
    {
      if ( d[j] < p )
      {
        k = j;
        p = d[j];
      }
    }

    if ( k != i )
    {
      d[k] = d[i];
      d[i] = p;

      for ( ii = 0; ii < n; ii++ )
      {
        t         = z[ii+i*n];
        z[ii+i*n] = z[ii+k*n];
        z[ii+k*n] = t;
      }
    }
  }

  return ierr;
}
/******************************************************************************/

int imtqlv ( int n, double d[], double e[], double e2[], double w[], int ind[] )

/******************************************************************************/
/*
  Purpose:

    IMTQLV computes all eigenvalues of a real symmetric tridiagonal matrix.

  Discussion:

    IMTQLV finds the eigenvalues of a symmetric tridiagonal matrix by 
    the implicit QL method and associates with them their corresponding 
    submatrix indices.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double D[N], the diagonal elements of the input matrix.

    Input, double E[N], the subdiagonal elements of the input matrix
    in E(2:N).  E(1) is arbitrary.

    Input/output, double E2[N].  On input, the squares of the
    corresponding elements of E.  E2(1) is arbitrary.  On output, elements of 
    E2 corresponding to elements of E regarded as negligible have been
    replaced by zero, causing the matrix to split into a direct sum of
    submatrices.  E2(1) is also set to zero.

    Output, double W[N], the eigenvalues in ascending order.  If an
    error break; is made, the eigenvalues are correct and ordered for
    indices 1,2,...IERR-1, but may not be the smallest eigenvalues.

    Output, int IND[N], the submatrix indices associated with 
    the corresponding eigenvalues in W: 1 for eigenvalues belonging to the
    first submatrix from the top, 2 for those belonging to the second
    submatrix, and so on.

    Output, int IMTQLV, error flag.
    0, for normal return,
    J, if the J-th eigenvalue has not been determined after 30 iterations.
*/
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ierr;
  int its;
  int k;
  int l;
  int m;
  double p;
  double r;
  double *rv1;
  double s;
  bool skip;
  int tag;
  double tst1;
  double tst2;

  ierr = 0;

  k = -1;
  tag = -1;
  for ( i = 0; i < n; i++ )
  {
    w[i] = d[i];
  }
  e2[0] = 0.0;

  rv1 = ( double * ) malloc ( n * sizeof ( double ) );
  for ( i = 0; i < n - 1; i++ )
  {
    rv1[i] = e[i+1];
  }
  rv1[n-1] = 0.0;

  for ( l = 0; l < n; l++ )
  {
    its = 0;
/*
  Look for a small sub-diagonal element.
*/
    while ( true )
    {
      for ( m = l; m < n; m++ )
      {
        if ( m == n - 1 )
        {
          break;
        }

        tst1 = fabs ( w[m] ) + fabs ( w[m+1] );
        tst2 = tst1 + fabs ( rv1[m] );

        if ( tst2 == tst1 )
        {
          break;
        }
/*
  Guard against underflowed element of E2.
*/
        if ( e2[m+1] == 0.0 )
        {
          k = m;
          tag = tag + 1;
          break;
        }
      }
 
      if ( k < m )
      {
        if ( m < n - 1 )
        {
          e2[m+1] = 0.0;
        }
        k = m;
        tag = tag + 1;
      }

      p = w[l];

      if ( m == l )
      {
        for ( i = l; 0 <= i; i-- )
        {
          if ( i == 0 )
          {
            w[i] = p;
            ind[i] = tag;
          }
          else if ( w[i-1] <= p )
          {
            w[i] = p;
            ind[i] = tag;
            break;
          }
          else
          {
            w[i] = w[i-1];
            ind[i] = ind[i-1];
          }
        }

        break;
      }
      else
      {
        if ( 30 <= its )
        {
          free ( rv1 );
          ierr = l + 1;
          return ierr;
        }

        its = its + 1;
/*
  Form shift.
*/
        g = ( w[l+1] - p ) / ( 2.0 * rv1[l] );
        r = pythag ( g, 1.0 );
        g = w[m] - p + rv1[l] / ( g + fabs ( r ) * r8_sign ( g ) );
        s = 1.0;
        c = 1.0;
        p = 0.0;

        skip = false;

        for ( i = m - 1; l <= i; i-- )
        {
          f = s * rv1[i];
          b = c * rv1[i];
          r = pythag ( f, g );
          rv1[i+1] = r;

          if ( r == 0.0 )
          {
            w[i+1] = w[i+1] - p;
            rv1[m] = 0.0;
            skip = true;
            break;
          }

          s = f / r;
          c = g / r;
          g = w[i+1] - p;
          r = ( w[i] - g ) * s + 2.0 * c * b;
          p = s * r;
          w[i+1] = g + p;
          g = c * r - b;
        }

        if ( ! skip )
        {
          w[l] = w[l] - p;
          rv1[l] = g;
          rv1[m] = 0.0;
        }
      }
    }
  }

  free ( rv1 );

  return ierr;
}
/******************************************************************************/

int invit ( int n, double a[], double wr[], double wi[], bool select[], 
  int mm, int *m, double z[] )

/******************************************************************************/
/*
  Purpose:

    INVIT is DUMMY CODE right now.
*/
{
  if ( true )
  {
    printf ( "\n" );
    printf ( "INVIT - Fatal error!\n" );
    printf ( "  This is just DUMMY CODE right now.\n" );
    exit ( 1 );
  }
  return 0;
}
/******************************************************************************/

int minfit ( int nm, int m, int n, double a[], double w[], int ip, double b[] )

/******************************************************************************/
/*
  Purpose:

    MINFIT: least squares problem for a real overdetermined linear system.

  Discussion:

    MINFIT is part of an algorithm for solving general linear
    systems of the form A*X=B.

    It determines the singular value decomposition
      A = U * S * V'
    of a real M by N rectangular matrix, forming U' * B
    rather than U.  Householder bidiagonalization and a variant of the
    QR algorithm are used.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 March 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int NM, the leading dimension of the
    two-dimensional arrays.  NM must be at least as large as the maximum
    of M and N.

    Input, int M, the number of rows of A and B.

    Input, int N, the number of columns of A, and the order of V.

    Input/output, double A[NM,N]. On input, the rectangular
    coefficient matrix.  On output, A has been overwritten by the orthogonal
    matrix V of the decomposition in its first N rows and columns.  If an
    error return is made, the columns of V corresponding to indices of correct
    singular values should be correct.

    Output, double W[N], the singular values of A.  These are the
    diagonal elements of S.  They are unordered.  If an error return is made, the
    singular values should be correct for indices IERR+1, IERR+2,...,N.

    Input, int IP, is the number of columns of B.  IP can be zero.

    Input/output, double B[NM,IP].  On input, the constant column
    matrix.  On output, B has been overwritten by U'*B.  If an error return is
    made, the rows of U'*B corresponding to indices of correct singular values
    should be correct.

    Output, int MINFIT, error flag.
    0, for normal return,
    K, if the K-th singular value has not been determined after 30 iterations.
*/
{
  double c;
  double f;
  double g;
  double h;
  int i;
  int ierr;
  int its;
  int j;
  int k;
  int l;
  double *rv1;
  double s;
  double scale;
  bool skip;
  double tst1;
  double tst2;
  double x;
  double y;
  double z;

  ierr = 0;

  rv1 = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Householder reduction to bidiagonal form.
*/
  g = 0.0;
  scale = 0.0;
  x = 0.0;

  for ( i = 0; i < n; i++ )
  {
    l = i + 1;
    rv1[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;

    if ( i < m )
    {
      scale = 0.0;
      for ( k = i; k < m; k++ )
      {
        scale = scale + ( fabs ( a[k+i*nm] ) );
      }

      if ( scale != 0.0 )
      {
        for ( k = i; k < m; k++ )
        {
          a[k+i*nm] = a[k+i*nm] / scale;
          s = s + pow ( a[k+i*nm], 2 );
        }
        f = a[i+i*nm];
        g = - sqrt ( s ) * r8_sign ( f );
        h = f * g - s;
        a[i+i*nm] = f - g;

        for ( j = i + 1; j < n; j++ )
        {
          s = 0.0;
          for ( k = i; k < m; k++ )
          {
            s = s + a[k+i*nm] * a[k+j*nm];
          }
          f = s / h;
          for ( k = i; k < m; k++ )
          {
            a[k+j*nm] = a[k+j*nm] + f * a[k+i*nm];
          }
        }

        for ( j = 0; j < ip; j++ )
        {
          s = 0.0;
          for ( k = i; k < m; k++ )
          {
            s = s + a[k+i*nm] * b[k+j*nm];
          }
          for ( k = i; k < m; k++ )
          {
            b[k+j*nm] = b[k+j*nm] + s * a[k+i*nm] / h;
          }
        }
        for ( k = i; k < m; k++ )
        {
          a[k+i*nm] = scale * a[k+i*nm];
        }
      }
    }

    w[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;

    if ( i < m && i != n - 1 )
    {
      for ( k = i + 1; k < n; k++ )
      {
        scale = scale + fabs ( a[i+k*nm] );
      }

      if ( scale != 0.0 )
      {
        for ( k = i + 1; k < n; k++ )
        {
          a[i+k*nm] = a[i+k*nm] / scale;
          s = s + pow ( a[i+k*nm], 2 );
        }
        f = a[i+(i+1)*nm];
        g = - sqrt ( s ) * r8_sign ( f );
        h = f * g - s;
        a[i+(i+1)*nm] = f - g;
        for ( k = i + 1; k < n; k++ )
        {
          rv1[k] = a[i+k*nm] / h;
        }

        for ( j = i + 1; j < m; j++ )
        {
          s = 0.0;
          for ( k = i + 1; k < n; k++ )
          {
            s = s + a[j+k*nm] * a[i+k*nm];
          }
          for ( k = i + 1; k < n; k++ )
          {
            a[j+k*nm] = a[j+k*nm] + s * rv1[k];
          }
        }
        for ( k = i + 1; k < n; k++ )
        {
          a[i+k*nm] = scale * a[i+k*nm];
        }
      }
    }
    x = r8_max ( x, fabs ( w[i] ) + fabs ( rv1[i] ) );
  }
/*
  Accumulation of right-hand transformations.
*/
  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( i < n - 1 )
    {
      if ( g != 0.0 )
      {
        for ( j = i + 1; j < n; j++ )
        {
          a[j+i*nm] = ( a[i+j*nm] / a[i+(i+1)*nm] ) / g;
        }

        for ( j = i + 1; j < n; j++ )
        {
          s = 0.0;
          for ( k = i + 1; k < n; k++ )
          {
            s = s + a[i+k*nm] * a[k+j*nm];
          }
          for ( k = i + 1; k < n; k++ )
          {
            a[k+j*nm] = a[k+j*nm] + s * a[k+i*nm];
          }
        }
      }
      for ( j = i + 1; j < n; j++ )
      {
        a[i+j*nm] = 0.0;
        a[j+i*nm] = 0.0;
      }
    }
    a[i+i*nm] = 1.0;
    g = rv1[i];
  }

  if ( m < n && ip != 0 )
  {
    for ( j = 0; j < ip; j++ )
    {
      for ( i = m + 1; i < n; i++ )
      {
        b[i+j*nm] = 0.0;
      }
    }
  }
/*
  Diagonalization of the bidiagonal form.
*/
  tst1 = x;

  for ( k = n - 1; 0 <= k; k-- )
  {
    its = 0;

    while ( true )
    {
/*
  Test for splitting.
*/
      skip = false;

      for ( l = k; 0 <= l; l-- )
      {
        tst2 = tst1 + fabs ( rv1[l] );

        if ( tst2 == tst1 )
        {
          skip = true;
          break;
        }

        tst2 = tst1 + fabs ( w[l-1] );

        if ( tst2 == tst1 )
        {
          break;
        }
      }
/*
  Cancellation of RV1[l] if l greater than 1.
*/
      if ( ! skip )
      {
        c = 0.0;
        s = 1.0;

        for ( i = l; i <= k; i++ )
        {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];
          tst2 = tst1 + fabs ( f );

          if ( tst2 == tst1 )
          {
            break;
          }

          g = w[i];
          h = pythag ( f, g );
          w[i] = h;
          c = g / h;
          s = - f / h;

          for ( j = 0; j < ip; j++ )
          {
            y = b[l-1+j*nm];
            z = b[i+j*nm];
            b[l-1+j*nm] =  y * c + z * s;
            b[i+j*nm] =  - y * s + z * c;
          }
        }
      }
/*
  Test for convergence.
*/
      z = w[k];

      if ( l == k )
      {
        if ( z < 0.0 )
        {
          w[k] = - z;
          for ( j = 0; j < n; j++ )
          {
            a[j+k*nm] = - a[j+k*nm];
          }
        }
        free ( rv1 );
        return ierr;
      }
/*
  Shift from bottom 2 by 2 minor.
*/
      if ( 30 <= its )
      {
        ierr = k;
        free ( rv1 );
        return ierr;
      }

      its = its + 1;
      x = w[l];
      y = w[k-1];
      g = rv1[k-1];
      h = rv1[k];
      f = 0.5 * ( ( ( g + z ) / h ) * ( ( g - z ) / y ) + y / h - h / y );
      g = pythag ( f, 1.0 );
      f = x - ( z / x ) * z + ( h / x ) * ( y / ( f + fabs ( g ) * r8_sign ( f ) ) - h );
/*
  Next QR transformation.
*/
      c = 1.0;
      s = 1.0;

      for ( i = l + 1; i <= k; i++ )
      {
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = pythag ( f, h );
        rv1[i-1] = z;
        c = f / z;
        s = h / z;
        f =   x * c + g * s;
        g = - x * s + g * c;
        h = y * s;
        y = y * c;

        for ( j = 0; j < n; j++ )
        {
          x = a[j+(i-1)*nm];
          z = a[j+i*nm];
          a[j+(i-1)*nm] =  x * c + z * s;
          a[j+i*nm] =  - x * s + z * c;
        }

        z = pythag ( f, h );
        w[i-1] = z;

        if ( z != 0.0 )
        {
          c = f / z;
          s = h / z;
        }

        f =   c * g + s * y;
        x = - s * g + c * y;

        for ( j = 0; j < ip; j++ )
        {
          y = b[i-1+j*nm];
          z = b[i+j*nm];
          b[i-1+j*nm] =  y * c + z * s;
          b[i+j*nm] =  - y * s + z * c;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }

  free ( rv1 );

  return ierr;
}

/******************************************************************************/

void ortbak ( int n, int low, int igh, double a[], double ort[], int m, 
  double z[] )

/******************************************************************************/
/*
  Purpose:

    ORTBAK determines eigenvectors by undoing the ORTHES transformation.

  Discussion:

    ORTBAK forms the eigenvectors of a real general
    matrix by back transforming those of the corresponding
    upper Hessenberg matrix determined by ORTHES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing
    routine BALANC.  If BALANC has not been used, set LOW = 1 and IGH equal
    to the order of the matrix.

    Input, double A[N,IGH], contains information about the 
    orthogonal transformations used in the reduction by ORTHES in its strict
    lower triangle.

    Input/output, double ORT[IGH], contains further information
    about the transformations used in the reduction by ORTHES.  On output, ORT
    has been altered.

    Input, int M, the number of columns of Z to be back
    transformed.

    Input/output, double Z[N,N].  On input, the real and imaginary
    parts of the eigenvectors to be back transformed in the first M columns.
    On output, the real and imaginary parts of the transformed eigenvectors.
*/
{
  double g;
  int i;
  int j;
  int mp;

  if ( m == 0 )
  {
    return;
  }

  for ( mp = igh - 1; low + 1 <= mp; mp-- )
  {
    if ( a[mp+(mp-1)*n] != 0.0 )
    {
      for ( i = mp + 1; i <= igh; i++ )
      {
        ort[i] = a[i+(mp-1)*n];
      }
      for ( j = 0; j < m; j++ )
      {
        g = 0.0;
        for ( i = mp; i <= igh; i++ )
        {
          g = g + ort[i] * z[i+j*n];
        }
        g = ( g / ort[mp] ) / a[mp+(mp-1)*n];
        for ( i = mp; i <= igh; i++ )
        {
          z[i+j*n] = z[i+j*n] + g * ort[i];
        }
      }
    }
  }

  return;
}
/******************************************************************************/

void orthes ( int n, int low, int igh, double a[], double ort[] )

/******************************************************************************/
/*
  Purpose:

    ORTHES transforms a real general matrix to upper Hessenberg form.

  Discussion:

    ORTHES is given a real general matrix, and reduces a submatrix
    situated in rows and columns LOW through IGH to upper Hessenberg form by
    orthogonal similarity transformations.

  Modified:

    03 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing
    routine BALANC.  If BALANC has not been used, set LOW = 1 and IGH = N.

    Input/output, double A[N,N].  On input, the matrix.  On output,
    the Hessenberg matrix.  Information about the orthogonal transformations
    used in the reduction is stored in the remaining triangle under the
    Hessenberg matrix.

    Output, double ORT[IGH], contains further information about the
    transformations.
*/
{
  double f;
  double g;
  double h;
  int i;
  int j;
  int m;
  double scale;

  for ( m = low + 1; m <= igh - 1; m++ )
  {
    h = 0.0;
    ort[m] = 0.0;
    scale = 0.0;
/*
  Scale the column.
*/
    for ( i = m; i <= igh; i++ )
    {
      scale = scale + fabs ( a[i+(m-1)*n] );
    }

    if ( scale != 0.0 )
    {
      for ( i = igh; m <= i; i-- )
      {
        ort[i] = a[i+(m-1)*n] / scale;
        h = h + ort[i] * ort[i];
      }

      g = - sqrt ( h ) * r8_sign ( ort[m] );
      h = h - ort[m] * g;
      ort[m] = ort[m] - g;
/*
  Form (I-(U*Ut)/h) * A.
*/
      for ( j = m; j < n; j++ )
      {
        f = 0.0;
        for ( i = igh; m <= i; i-- )
        {
          f = f + ort[i] * a[i+j*n];
        }
        f = f / h;

        for ( i = m; i <= igh; i++ )
        {
          a[i+j*n] = a[i+j*n] - f * ort[i];
        }
      }
/*
  Form (I-(u*ut)/h) * A * (I-(u*ut)/h).
*/
      for ( i = 0; i <= igh; i++ )
      {
        f = 0.0;
        for ( j = igh; m <= j; j-- )
        {
          f = f + ort[j] * a[i+j*n];
        }
        for ( j = m; j <= igh; j++ )
        { 
          a[i+j*n] = a[i+j*n] - f * ort[j] / h;
        }
      }
      ort[m] = scale * ort[m];
      a[m+(m-1)*n] = scale * g;
    }
  }

  return;
}
/******************************************************************************/

void ortran ( int n, int low, int igh, double a[], double ort[], double z[] )

/******************************************************************************/
/*
  Purpose:

    ORTRAN accumulates similarity transformations generated by ORTHES.

  Discussion:

    ORTRAN accumulates the orthogonal similarity
    transformations used in the reduction of a real general
    matrix to upper Hessenberg form by ORTHES.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int LOW, IGH, are determined by the balancing
    routine BALANC.  If BALANC has not been used, set LOW = 1, IGH = N.

    Input, double A[N,IGH], contains information about the 
    orthogonal transformations used in the reduction by ORTHES in its strict 
    lower triangle.

    Input/output, double ORT[IGH], contains further information
    about the transformations used in the reduction by ORTHES.  On output, ORT
    has been further altered.

    Output, double Z[N,N], contains the transformation matrix
    produced in the reduction by ORTHES.
*/
{
  double g;
  int i;
  int j;
  int mp;
/*
  Initialize Z to the identity matrix.
*/
  r8mat_identity ( n, z );

  if ( igh - low < 2 )
  {
    return;
  }

  for ( mp = igh - 1; low + 1 <= mp; mp-- )
  {
    if ( a[mp+(mp-1)*n] != 0.0 )
    {
      for ( i = mp + 1; i <= igh; i++ )
      {
        ort[i] = a[i+(mp-1)*n];
      }
      for ( j = mp; j <= igh; j++ )
      {
        g = 0.0;
        for ( i = mp; i <= igh; i++ )
        {
          g = g + ort[i] * z[i+j*n];
        }
        g = ( g / ort[mp] ) / a[mp+(mp-1)*n];

        for ( i = mp; i <= igh; i++ )
        {
          z[i+j*n] = z[i+j*n] + g * ort[i];
        }
      }
    }
  }

  return;
}
/******************************************************************************/

double pythag ( double a, double b )

/******************************************************************************/
/*
  Purpose:

    PYTHAG computes SQRT ( A * A + B * B ) carefully.

  Discussion:

    The formula

      PYTHAG = sqrt ( A * A + B * B )

    is reasonably accurate, but can fail if, for example, A^2 is larger
    than the machine overflow.  The formula can lose most of its accuracy
    if the sum of the squares is very large or very small.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Modified:

    08 November 2012

  Parameters:

    Input, double A, B, the two legs of a right triangle.

    Output, double PYTHAG, the length of the hypotenuse.
*/
{
  double p;
  double r;
  double s;
  double t;
  double u;

  p = r8_max ( fabs ( a ), fabs ( b ) );

  if ( p != 0.0 )
  {
    r = r8_min ( fabs ( a ), fabs ( b ) ) / p;
    r = r * r;

    while ( 1 )
    {
      t = 4.0 + r;

      if ( t == 4.0 )
      {
        break;
      }

      s = r / t;
      u = 1.0 + 2.0 * s;
      p = u * p;
      r = ( s / u ) * ( s / u ) * r;
    }
  }
  return p;
}
/******************************************************************************/

void qzhes ( int n, double a[], double b[], bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    QZHES carries out transformations for a generalized eigenvalue problem.

  Discussion:

    QZHES is the first step of the QZ algorithm
    for solving generalized matrix eigenvalue problems.

    QZHES accepts a pair of real general matrices and
    reduces one of them to upper Hessenberg form and the other
    to upper triangular form using orthogonal transformations.
    it is usually followed by QZIT, QZVAL and, possibly, QZVEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrices.

    Input/output, double A[N,N].  On input, the first real general
    matrix.  On output, A has been reduced to upper Hessenberg form.  The
    elements below the first subdiagonal have been set to zero.

    Input/output, double B[N,N].  On input, a real general matrix.
    On output, B has been reduced to upper triangular form.  The elements
    below the main diagonal have been set to zero.

    Input, bool MATZ, is true if the right hand transformations
    are to be accumulated for later use in computing eigenvectors.

    Output, double Z[N,N], contains the product of the right hand
    transformations if MATZ is true.
*/
{
  int i;
  int j;
  int k;
  int l;
  double r;
  double rho;
  double s;
  double t;
  double u1;
  double u2;
  double v1;
  double v2;
/*
  Set Z to the identity matrix.
*/
  if ( matz )
  {
    r8mat_identity ( n, z );
  }
/*
  Reduce B to upper triangular form.
*/
  if ( n <= 1 )
  {
    return;
  }

  for ( l = 0; l < n - 1; l++ )
  {
    s = 0.0;
    for ( i = l + 1; i < n; i++ )
    {
      s = s + fabs ( b[i+l*n] );
    }

    if ( s != 0.0 )
    {
      s = s + fabs ( b[l+l*n] );
      for ( i = l; i < n; i++ )
      {
        b[i+l*n] = b[i+l*n] / s;
      }
      r = 0.0;
      for ( i = l; i < n; i++ )
      {
        r = r + b[i+l*n] * b[i+l*n];
      }
      r = sqrt ( r );
      r = fabs ( r ) * r8_sign ( b[l+l*n] );
      b[l+l*n] = b[l+l*n] + r;
      rho = r * b[l+l*n];

      for ( j = l + 1; j < n; j++ )
      {
        t = 0.0;
        for ( i = l; i < n; i++ )
        {
          t = t + b[i+l*n] * b[i+j*n];
        }
        for ( i = l; i < n; i++ )
        {
          b[i+j*n] = b[i+j*n] - t * b[i+l*n] / rho;
        }
      }

      for ( j = 0; j < n; j++ )
      {
        t = 0.0;
        for ( i = l; i < n; i++ )
        {
          t = t + b[i+l*n] * a[i+j*n];
        }
        for ( i = l; i < n; i++ )
        {
          a[i+j*n] = a[i+j*n] - t * b[i+l*n] / rho;
        }
      }

      b[l+l*n] = - s * r;
      for ( i = l + 1; i < n; i++ )
      {
        b[i+l*n] = 0.0;
      }
    }
  }
/*
  Reduce A to upper Hessenberg form, while keeping B triangular.
*/
  for ( k = 0; k < n - 2; k++ )
  {
    for ( l = n - 2; k + 1 <= l; l-- )
    {
/*
  Zero A[l+1+k*n].
*/
      s = fabs ( a[l+k*n] ) + fabs ( a[l+1+k*n] );

      if ( s != 0.0 )
      {
        u1 = a[l+k*n] / s;
        u2 = a[l+1+k*n] / s;
        r = sqrt ( u1 * u1 + u2 * u2 ) * r8_sign ( u1 );
        v1 = - ( u1 + r ) / r;
        v2 = - u2 / r;
        u2 = v2 / v1;

        for ( j = k; j < n; j++ )
        {
          t = a[l+j*n] + u2 * a[l+1+j*n];
          a[l+j*n] = a[l+j*n] + t * v1;
          a[l+1+j*n] = a[l+1+j*n] + t * v2;
        }

        a[l+1+k*n] = 0.0;

        for ( j = l; j < n; j++ )
        {
          t = b[l+j*n] + u2 * b[l+1+j*n];
          b[l+j*n] = b[l+j*n] + t * v1;
          b[l+1+j*n] = b[l+1+j*n] + t * v2;
        }
/*
  Zero B[l+1+l*n].
*/
        s = fabs ( b[l+1+(l+1)*n] ) + fabs ( b[l+1+l*n] );

        if ( s != 0.0 )
        {
          u1 = b[l+1+(l+1)*n] / s;
          u2 = b[l+1+l*n] / s;
          r = sqrt ( u1 * u1 + u2 * u2 ) * r8_sign ( u1 );
          v1 =  - ( u1 + r ) / r;
          v2 = - u2 / r;
          u2 = v2 / v1;

          for ( i = 0; i <= l + 1; i++ )
          {
            t = b[i+(l+1)*n] + u2 * b[i+l*n];
            b[i+(l+1)*n] = b[i+(l+1)*n] + t * v1;
            b[i+l*n] = b[i+l*n] + t * v2;
          }

          b[l+1+l*n] = 0.0;

          for ( i = 0; i < n; i++ )
          {
            t = a[i+(l+1)*n] + u2 * a[i+l*n];
            a[i+(l+1)*n] = a[i+(l+1)*n] + t * v1;
            a[i+l*n] = a[i+l*n] + t * v2;
          }

          if ( matz )
          {
            for ( i = 0; i < n; i++ )
            {
              t = z[i+(l+1)*n] + u2 * z[i+l*n];
              z[i+(l+1)*n] = z[i+(l+1)*n] + t * v1;
              z[i+l*n] = z[i+l*n] + t * v2;
            }
          }
        }
      }
    }
  }

  return;
}
/******************************************************************************/

int qzit ( int n, double a[], double b[], double eps1, bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    QZIT is DUMMY CODE right now.
*/
{
  if ( true )
  {
    printf ( "\n" );
    printf ( "QZIT - Fatal error!\n" );
    printf ( "  This is just DUMMY CODE right now.\n" );
    exit ( 1 );
  }
  return 0;
}
/******************************************************************************/

void qzval ( int n, double a[], double b[], double alfr[], double alfi[], 
  double beta[], bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    QZVAL computes eigenvalues for a generalized eigenvalue problem.

  Discussion:

    QZVAL is the third step of the QZ algorithm
    for solving generalized matrix eigenvalue problems.

    QZVAL accepts a pair of real matrices, one of them
    in quasi-triangular form and the other in upper triangular form.
    It reduces the quasi-triangular matrix further, so that any
    remaining 2-by-2 blocks correspond to pairs of complex
    eigenvalues, and returns quantities whose ratios give the
    generalized eigenvalues.  It is usually preceded by QZHES
    and QZIT and may be followed by QZVEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrices.

    Input/output, double A[N*N].  On input, a real upper
    quasi-triangular matrix.  On output, A has been reduced further to a
    quasi-triangular matrix in which all nonzero subdiagonal elements
    correspond to pairs of complex eigenvalues.

    Input/output, double B[N*N].  On input, a real upper triangular 
    matrix.  In addition, location B(n,1) contains the tolerance quantity EPSB
    computed and saved in QZIT.  On output, B is still in upper triangular
    form, although its elements have been altered.  B(N,1) is unaltered.

    Output, double ALFR[N], ALFI[N], the real and imaginary parts of
    the diagonal elements of the triangular matrix that would be obtained
    if A were reduced completely to triangular form by unitary
    transformations.  Non-zero values of ALFI occur in pairs, the first
    member positive and the second negative.

    Output, double BETA[N], the diagonal elements of the 
    corresponding B, normalized to be real and non-negative.  The generalized 
    eigenvalues are then the ratios (ALFR + I * ALFI) / BETA.

    Input, bool MATZ, should be true if the right hand transformations
    are to be accumulated for later use in computing eigenvectors, and
    false otherwise.

    Input/output, double Z[N*N], is only used if MATZ is true.
    On input, the transformation matrix produced in the reductions by QZHES
    and QZIT, if performed, or else the identity matrix.  On output,
    the product of the right hand transformations for all three steps.
*/
{
  double a1;
  double a11;
  double a11i;
  double a11r;
  double a12;
  double a12i;
  double a12r;
  double a1i;
  double a2;
  double a21;
  double a22;
  double a22i;
  double a22r;
  double a2i;
  double an = 0.0 ;
  double b11;
  double b12;
  double b22;
  double bn;
  double c;
  double cq;
  double cz;
  double d;
  double di;
  double dr;
  double e = 0.0;
  double ei;
  int en;
  double epsb;
  int i;
  int isw;
  int j;
  int na;
  double r;
  double s;
  double sqi;
  double sqr;
  double ssi;
  double ssr;
  double szi;
  double szr;
  double t;
  double ti;
  double tr;
  double u1;
  double u2;
  double v1;
  double v2;

  epsb = b[n-1+0*n];
  isw = 1;
/*
  Find eigenvalues of quasi-triangular matrices.
*/
  for ( en = n - 1; 0 <= en; en-- )
  {
    na = en - 1;

    if ( isw == 2 )
    {
      isw = 3 - isw;
      continue;
    }
/*
  1-by-1 block, one real root.
*/
    if ( en == 0 || a[en+na*n] == 0.0 )
    {
      alfr[en] = a[en+en*n];
      if ( b[en+en*n] < 0.0 )
      {
        alfr[en] = - alfr[en];
      }
      beta[en] = fabs ( b[en+en*n] );
      alfi[en] = 0.0;
      continue;
    }
/*
  2-by-2 block.
*/
    if ( fabs ( b[na+na*n] ) <= epsb )
    {
      a1 = a[na+na*n];
      a2 = a[en+na*n];
    }
    else
    {
      if ( fabs ( b[en+en*n] ) <= epsb )
      {
        a1 = a[en+en*n];
        a2 = a[en+na*n];
        bn = 0.0;
      }
      else
      {
        an = fabs ( a[na+na*n] ) + fabs ( a[na+en*n] ) + fabs ( a[en+na*n] ) 
          + fabs ( a[en+en*n] );
        bn = fabs ( b[na+na*n] ) + fabs ( b[na+en*n] ) + fabs ( b[en+en*n] );
        a11 = a[na+na*n] / an;
        a12 = a[na+en*n] / an;
        a21 = a[en+na*n] / an;
        a22 = a[en+en*n] / an;
        b11 = b[na+na*n] / bn;
        b12 = b[na+en*n] / bn;
        b22 = b[en+en*n] / bn;
        e = a11 / b11;
        ei = a22 / b22;
        s = a21 / ( b11 * b22 );
        t = ( a22 - e * b22 ) / b22;

        if ( fabs ( ei ) < fabs ( e ) )
        {
          e = ei;
          t = ( a11 - e * b11 ) / b11;
        }

        c = 0.5 * ( t - s * b12 );
        d = c * c + s * ( a12 - e * b12 );
/*
  Two complex roots.
*/
        if ( d < 0.0 )
        {
          e = e + c;
          ei = sqrt ( - d );
          a11r = a11 - e * b11;
          a11i = ei * b11;
          a12r = a12 - e * b12;
          a12i = ei * b12;
          a22r = a22 - e * b22;
          a22i = ei * b22;

          if ( fabs ( a21 ) + fabs ( a22r ) + fabs ( a22i ) <= 
               fabs ( a11r ) + fabs ( a11i ) + fabs ( a12r ) + fabs ( a12i ) )
          {
            a1 = a12r;
            a1i = a12i;
            a2 = - a11r;
            a2i = - a11i;
          }
          else
          {
            a1 = a22r;
            a1i = a22i;
            a2 = - a21;
            a2i = 0.0;
          }
/*
  Choose complex Z.
*/
          cz = sqrt ( a1 * a1 + a1i * a1i );

          if ( cz != 0.0 )
          {
            szr = ( a1 * a2 + a1i * a2i ) / cz;
            szi = ( a1 * a2i - a1i * a2 ) / cz;
            r = sqrt ( cz * cz + szr * szr + szi * szi );
            cz = cz / r;
            szr = szr / r;
            szi = szi / r;
          }
          else
          {
            szr = 1.0;
            szi = 0.0;
          }

          if ( ( fabs ( e ) + ei ) * bn <= an )
          {
            a1 = cz * b11 + szr * b12;
            a1i = szi * b12;
            a2 = szr * b22;
            a2i = szi * b22;
          }
          else
          {
            a1 = cz * a11 + szr * a12;
            a1i = szi * a12;
            a2 = cz * a21 + szr * a22;
            a2i = szi * a22;
          }
/*
  Choose complex Q.
*/
          cq = sqrt ( a1 * a1 + a1i * a1i );

          if ( cq != 0.0 )
          {
            sqr = ( a1 * a2 + a1i * a2i ) / cq;
            sqi = ( a1 * a2i - a1i * a2 ) / cq;
            r = sqrt ( cq * cq + sqr * sqr + sqi * sqi );
            cq = cq / r;
            sqr = sqr / r;
            sqi = sqi / r;
          }
          else
          {
            sqr = 1.0;
            sqi = 0.0;
          }
/*
  Compute diagonal elements that would result if transformations were applied.
*/
          ssr = sqr * szr + sqi * szi;
          ssi = sqr * szi - sqi * szr;
          i = 1;
          tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22;
          ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22;
          dr = cq * cz * b11 + cq * szr * b12 + ssr * b22;
          di = cq * szi * b12 + ssi * b22;

          while ( true )
          {
            t = ti * dr - tr * di;
 
            if ( t < 0.0 )
            {
              j = en;
            }
            else
            {
              j = na;
            }

            r = sqrt ( dr * dr + di * di );
            beta[j] = bn * r;
            alfr[j] = an * ( tr * dr + ti * di ) / r;
            alfi[j] = an * t / r;

            if ( i != 1 )
            {
              break;
            }

            i = 2;
            tr =   ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22;
            ti = - ssi * a11 - sqi * cz * a12 + cq * szi * a21;
            dr =   ssr * b11 - sqr * cz * b12 + cq * cz * b22;
            di = - ssi * b11 - sqi * cz * b12;
          }

          isw = 3 - isw;
          continue;
        }
/*
  Two real roots.
  Zero both A(EN,NA) and B(EN,NA).
*/
        e = e + ( c + sqrt ( d ) * r8_sign ( c ) );
        a11 = a11 - e * b11;
        a12 = a12 - e * b12;
        a22 = a22 - e * b22;

        if ( fabs ( a21 ) + fabs ( a22 ) <= fabs ( a11 ) + fabs ( a12 ) )
        {
          a1 = a12;
          a2 = a11;
        }
        else
        {
          a1 = a22;
          a2 = a21;
        }
      }
/*
  Choose and apply real Z.
*/
      s = fabs ( a1 ) + fabs ( a2 );
      u1 = a1 / s;
      u2 = a2 / s;
      r = sqrt ( u1 * u1 + u2 * u2 ) * r8_sign ( u1 );
      v1 = - ( u1 + r ) / r;
      v2 = - u2 / r;
      u2 = v2 / v1;

      for ( i = 0; i < en; i++ )
      {
        t = a[i+en*n] + u2 * a[i+na*n];
        a[i+en*n] = a[i+en*n] + t * v1;
        a[i+na*n] = a[i+na*n] + t * v2;
        t = b[i+en*n] + u2 * b[i+na*n];
        b[i+en*n] = b[i+en*n] + t * v1;
        b[i+na*n] = b[i+na*n] + t * v2;
      }

      if ( matz )
      {
        for ( i = 0; i < n; i++ )
        {
          t = z[i+en*n] + u2 * z[i+na*n];
          z[i+en*n] = z[i+en*n] + t * v1;
          z[i+na*n] = z[i+na*n] + t * v2;
        }
      }

      if ( bn == 0.0 )
      {
        a[en+na*n] = 0.0;
        b[en+na*n] = 0.0;
        alfr[na] = a[na+na*n];
        alfr[en] = a[en+en*n];
        if ( b[na+na*n] < 0.0 )
        {
          alfr[na] = - alfr[na];
        }

        if ( b[en+en*n] < 0.0 )
        {
          alfr[en] = - alfr[en];
        }

        beta[na] = fabs ( b[na+na*n] );
        beta[en] = fabs ( b[en+en*n] );
        alfi[en] = 0.0;
        alfi[na] = 0.0;
        isw = 3 - isw;
        continue;
      }

      if ( fabs ( e ) * bn <= an )
      {
        a1 = b[na+na*n];
        a2 = b[en+na*n];
      }
      else
      {
        a1 = a[na+na*n];
        a2 = a[en+na*n];
      }
    }
/*
  Choose and apply real Q.
*/
    s = fabs ( a1 ) + fabs ( a2 );

    if ( s != 0.0 )
    {
      u1 = a1 / s;
      u2 = a2 / s;
      r = sqrt ( u1 * u1 + u2 * u2 ) * r8_sign ( u1 );
      v1 = - ( u1 + r ) / r;
      v2 = - u2 / r;
      u2 = v2 / v1;

      for ( j = na; j < n; j++ )
      {
        t = a[na+j*n] + u2 * a[en+j*n];
        a[na+j*n] = a[na+j*n] + t * v1;
        a[en+j*n] = a[en+j*n] + t * v2;
        t = b[na+j*n] + u2 * b[en+j*n];
        b[na+j*n] = b[na+j*n] + t * v1;
        b[en+j*n] = b[en+j*n] + t * v2;
      }
    }

    a[en+na*n] = 0.0;
    b[en+na*n] = 0.0;
    alfr[na] = a[na+na*n];
    alfr[en] = a[en+en*n];
    if ( b[na+na*n] < 0.0 ) 
    {
      alfr[na] = - alfr[na];
    }

    if ( b[en+en*n] < 0.0 )
    {
      alfr[en] = - alfr[en];
    }

    beta[na] = fabs ( b[na+na*n] );
    beta[en] = fabs ( b[en+en*n] );
    alfi[en] = 0.0;
    alfi[na] = 0.0;
    isw = 3 - isw;
  }

  b[n-1+0*n] = epsb;

  return;
}
/******************************************************************************/

void qzvec ( int n, double a[], double b[], double alfr[], double alfi[], 
  double beta[], double z[] )

/******************************************************************************/
/*
  Purpose:

    QZVEC is DUMMY CODE right now.
*/
{
  if ( true )
  {
    printf ( "\n" );
    printf ( "QZVEC - Fatal error!\n" );
    printf ( "  This is just DUMMY CODE right now.\n" );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

double r8_epsilon ( )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON returns the R8 round off unit.

  Discussion:

    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON, the R8 round-off unit.
*/
{
  const double value = 2.220446049250313E-016;

  return value;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
/******************************************************************************/

double r8_sign ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN returns the sign of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose sign is desired.

    Output, double R8_SIGN, the sign of X.
*/
{
  double value;

  if ( x < 0.0 )
  {
    value = - 1.0;
  }
  else
  {
    value = + 1.0;
  }
  return value;
}
/******************************************************************************/

void r8mat_copy ( int m, int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_COPY copies one R8MAT to another.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A1[M*N], the matrix to be copied.

    Output, double A2[M*N], the copy of A1.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }

  return;
}
/******************************************************************************/

void r8mat_identity  ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_IDENTITY sets an R8MAT to the identity matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Output, double A[N*N], the N by N identity matrix.
*/
{
  int i;
  int j;
  int k;

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MM_NEW multiplies two matrices.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.

    Output, double R8MAT_MM[N1*N3], the product matrix C = A * B.
*/
{
  double *c;
  int i;
  int j;
  int k;

  c = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
/******************************************************************************/

double *r8mat_mmt_new ( int n1, int n2, int n3, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MMT_NEW computes C = A * B'.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N3*N2], the matrices to multiply.

    Output, double R8MAT_MMT_NEW[N1*N3], the product matrix C = A * B'.
*/
{
  double *c;
  int i;
  int j;
  int k;

  c = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[j+k*n3];
      }
    }
  }

  return c;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8's, which
    may be stored as a vector in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

double *r8mat_uniform_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_01_NEW fills an R8MAT with pseudorandom values scaled to [0,1].

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2009

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Philip Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0, otherwise the output value of SEED
    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
    been updated.

    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int j;
  int k;
  double *r;

  r = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

void r8mat_zeros ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_ZEROS zeroes an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, double A[M*N], a matrix of zeroes.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return;
}
/******************************************************************************/

void r8vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

void r8vec2_print ( int n, double a1[], double a2[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_PRINT prints an R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A1[N], double A2[N], the vectors to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %4d: %14g  %14g\n", i, a1[i], a2[i] );
  }

  return;
}
/******************************************************************************/

int ratqr ( int n, double eps1, double d[], double e[], double e2[], int m, 
  double w[], int ind[], double bd[], bool type, int idef )

/******************************************************************************/
/*
  Purpose:

    RATQR computes selected eigenvalues of a real symmetric tridiagonal matrix.

  Discussion:

    RATQR finds the algebraically smallest or largest eigenvalues of a 
    symmetric tridiagonal matrix by the rational QR method with Newton 
    corrections.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double EPS1.  On input, a theoretical absolute
    error tolerance for the computed eigenvalues.  If the input EPS1 is
    non-positive, or indeed smaller than its default value, it is reset at
    each iteration to the respective default value, namely, the product of
    the relative machine precision and the magnitude of the current eigenvalue
    iterate.  The theoretical absolute error in the K-th eigenvalue is usually
    not greater than K times EPS1.  On output, EPS1 is unaltered unless it has
    been reset to its (last) default value.

    Input, double D[N], the diagonal elements of the input matrix.

    Input, double E[N], the subdiagonal elements of the input matrix
    in E(2:N).  E(1) is arbitrary.

    Input/output, double E2[N].  On input, E2(2:N-1) contains the
    squares of the corresponding elements of E, and E2(1) is arbitrary.  On
    output, elements of E2 corresponding to elements of E regarded as
    negligible have been replaced by zero, causing the matrix to split into
    a direct sum of submatrices.  E2(1) is set to 0.0D+00 if the smallest
    eigenvalues have been found, and to 2.0D+00 if the largest eigenvalues
    have been found.  E2 is otherwise unaltered (unless overwritten by BD).

    Input, int M, the number of eigenvalues to be found.

    Output, double W[M], the M algebraically smallest eigenvalues in
    ascending order, or the M largest eigenvalues in descending order.
    If an error exit is made because of an incorrect specification of IDEF,
    no eigenvalues are found.  If the Newton iterates for a particular
    eigenvalue are not monotone, the best estimate obtained is returned
    and IERR is set.  W may coincide with D.

    Output, int IND[N], contains in its first M positions the submatrix
    indices associated with the corresponding eigenvalues in W:
    1 for eigenvalues belonging to the first submatrix from the top, 2 for
    those belonging to the second submatrix, and so on.

    Output, double BD[N], contains refined bounds for the
    theoretical errors of the corresponding eigenvalues in W.  These bounds
    are usually within the tolerance specified by EPS1.  BD may coincide
    with E2.

    Input, bool TYPE, should be set to TRUE if the smallest eigenvalues
    are to be found, and to FALSE if the largest eigenvalues are to be found.

    Input, int IDEF, should be set to 1 if the input matrix
    is known to be positive definite, to -1 if the input matrix is known to
    be negative  definite, and to 0 otherwise.

    Output, int RATQR, error flag.
    0, for normal return,
    6*N+1, if IDEF is set to 1 and TYPE to .true. when the matrix is not
      positive definite, or if IDEF is set to -1 and TYPE to .false.
      when the matrix is not negative definite,
    5*N+K, if successive iterates to the K-th eigenvalue are not monotone
      increasing, where K refers to the last such occurrence.
*/
{
  double delta;
  double ep;
  double err;
  double f;
  int i;
  int ierr;
  int ii;
  bool irreg;
  int j;
  int jdef;
  int k;
  double p;
  double q;
  double qp = 0.0 ;
  double r;
  double s;
  double tot;

  ierr = 0;
  jdef = idef;
  for ( i = 0; i < n; i++ )
  {
    w[i] = d[i];
  }

  if ( ! type )
  {
    for ( i = 0; i < n; i++ )
    {
      w[i] = - w[i];
    }
    jdef = - jdef;
  }

  while ( true )
  {
    err = 0.0;
    s = 0.0;
/*
  Look for small sub-diagonal entries and define initial shift
  from lower Gerschgorin bound.

  Copy E2 array into BD.
*/
    tot = w[0];
    q = 0.0;
    j = -1;

    for ( i = 0; i < n; i++ )
    {
      p = q;

      if ( i == 0 )
      {
        e2[i] = 0.0;
      }
      else if ( p <= ( fabs ( d[i] ) + fabs (  d[i-1] ) ) * r8_epsilon ( ) )
      {
        e2[i] = 0.0;
      }

      bd[i] = e2[i];
/*
  Count also if element of E2 has underflowed.
*/
      if ( e2[i] == 0.0 )
      {
        j = j + 1;
      }

      ind[i] = j;
      q = 0.0;
      if ( i < n - 1 )
      {
        q = fabs ( e[i+1] );
      }

      tot = r8_min ( tot, w[i] - p - q );
    }

    if ( jdef == 1 && tot < 0.0 )
    {
      tot = 0.0;
    }
    else
    {
      for ( i = 0; i < n; i++ )
      {
        w[i] = w[i] - tot;
      }
    }

    for ( k = 0; k < m; k++ )
    {
/*
  Next QR transformation.
*/
      irreg = true;

      while ( true )
      {
        tot = tot + s;
        delta = w[n-1] - s;
        i = n - 1;
        f = fabs ( tot ) * r8_epsilon ( );
        eps1 = r8_max ( eps1, f );

        if ( delta <= eps1 )
        {
          if ( delta < - eps1 )
          {
            ierr = 6 * n + 1;
            return ierr;
          }

          irreg = false;
          break;
        }
/*
  Replace small sub-diagonal squares by zero to reduce the incidence of
  underflows.
*/
        for ( j = k + 1; j < n; j++ )
        {
          if ( bd[j] <= r8_epsilon ( ) * r8_epsilon ( ) )
          {
            bd[j] = 0.0;
          }
        }

        f = bd[n-1] / delta;
        qp = delta + f;
        p = 1.0;

        for ( i = n - 2; k <= i; i-- )
        {
          q = w[i] - s - f;
          r = q / qp;
          p = p * r + 1.0;
          ep = f * r;
          w[i+1] = qp + ep;
          delta = q - ep;

          if ( delta <= eps1 )
          {
            if ( delta < - eps1 )
            {
              ierr = 6 * n + 1;
              return ierr;
            }
            irreg = false;
            break;
          }

          f = bd[i] / q;
          qp = delta + f;
          bd[i+1] = qp * ep;
        }

        if ( ! irreg )
        {
          break;
        }

        w[k] = qp;
        s = qp / p;

        if ( tot + s <= tot )
        {
          break;
        }
      }
/*
  Set error: irregular end of iteration.
  Deflate minimum diagonal element.
*/
      if ( irreg )
      {
        ierr = 5 * n + k;
        s = 0.0;
        delta = qp;

        for ( j = k; j < n; j++ )
        {
          if ( w[j] <= delta )
          {
            i = j;
            delta = w[j];
          }
        }
      }
/*
  Convergence.
*/
      if ( i < n - 1 )
      {
        bd[i+1] = bd[i] * f / qp;
      }

      ii = ind[i];

      for ( j = i - 1; k <= j; j-- )
      {
        w[j+1] = w[j] - s;
        bd[j+1] = bd[j];
        ind[j+1] = ind[j];
      }

      w[k] = tot;
      err = err + fabs ( delta );
      bd[k] = err;
      ind[k] = ii;
    }

    if ( type )
    {
      return ierr;
    }

    f = bd[0];
    e2[0] = 2.0;
    bd[0] = f;
/*
  Negate elements of W for largest values.
*/
    for ( i = 0; i < n; i++ )
    {
      w[i] = - w[i];
    }
    jdef = - jdef;

    break;
  }

  return ierr;
}
/******************************************************************************/

void rebak ( int n, double b[], double dl[], int m, double z[] )

/******************************************************************************/
/*
  Purpose:

    REBAK determines eigenvectors by undoing the REDUC transformation.

  Discussion:

    REBAK forms the eigenvectors of a generalized
    symmetric eigensystem by back transforming those of the
    derived symmetric matrix determined by REDUC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double B[N,N], contains information about the similarity
    transformation (Cholesky decomposition) used in the reduction by REDUC
    in its strict lower triangle.

    Input, double DL[N], further information about the 
    transformation.

    Input, int M, the number of eigenvectors to be back
    transformed.

    Input/output, double Z[N,M].  On input, the eigenvectors to be
    back transformed in its first M columns.  On output, the transformed
    eigenvectors.
*/
{
  double dot;
  int i;
  int j;
  int k;

  for ( j = 0; j < m; j++ )
  {
    for ( i = n - 1; 0 <= i; i-- )
    {
      dot = 0.0;
      for ( k = i + 1; k < n; k++ )
      {
        dot = dot + b[k+i*n] * z[k+j*n];
      }
      z[i+j*n] = ( z[i+j*n] - dot ) / dl[i];
    }
  }

  return;
}
/******************************************************************************/

void rebakb ( int n, double b[], double dl[], int m, double z[] )

/******************************************************************************/
/*
  Purpose:

    REBAKB determines eigenvectors by undoing the REDUC2 transformation.

  Discussion:

    REBAKB forms the eigenvectors of a generalized symmetric eigensystem by 
    back transforming those of the derived symmetric matrix determined 
    by REDUC2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double B[N,N], contains information about the similarity
    transformation (Cholesky decomposition) used in the reduction by REDUC2
    in its strict lower triangle.

    Input, double DL[N], further information about the 
    transformation.

    Input, int M, the number of eigenvectors to be back
    transformed.

    Input/output, double Z[N,M].  On input, the eigenvectors to be 
    back transformed in its first M columns.  On output, the transformed
    eigenvectors.
*/
{
  int i;
  int j;
  int k;
  double t;

  for ( j = 0; j < m; j++ )
  {
    for ( i = n - 1; 0 <= i; i-- )
    {
      t = dl[i] * z[i+j*n];
      for ( k = 0; k < i; k++ )
      {
        t = t + b[i+k*n] * z[k+j*n];
      }
      z[i+j*n] = t;
    }
  }

  return;
}
/******************************************************************************/

int reduc ( int n, double a[], double b[], double dl[] )

/******************************************************************************/
/*
  Purpose:

    REDUC reduces the eigenvalue problem A*x=lambda*B*x to A*x=lambda*x.

  Discussion:

    REDUC reduces the generalized symmetric eigenproblem
    a x=(lambda) b x, where B is positive definite, to the standard
    symmetric eigenproblem using the Cholesky factorization of B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrices A and B.  If the
    Cholesky factor L of B is already available, N should be prefixed with a
    minus sign.

    Input/output, double A[N,N].  On input, A contains a real
    symmetric matrix.  Only the full upper triangle of the matrix need be
    supplied.  On output, A contains in its full lower triangle the full lower
    triangle of the symmetric matrix derived from the reduction to the
    standard form.  The strict upper triangle of a is unaltered.

    Input/output, double B[N,N].  On input, the real symmetric
    input matrix.  Only the full upper triangle of the matrix need be supplied.
    If N is negative, the strict lower triangle of B contains, instead, the
    strict lower triangle of its Cholesky factor L.  In any case, on output,
    B contains in its strict lower triangle the strict lower triangle of
    its Cholesky factor L.  The full upper triangle of B is unaltered.

    Input/output, double DL[N].  If N is negative, then the DL
    contains the diagonal elements of L on input.  In any case, DL will contain
    the diagonal elements of L on output,

    Output, int REDUC, error flag.
    0, for normal return,
    7*N+1, if B is not positive definite.
*/
{
  int i;
  int ierr;
  int j;
  int k;
  int nn;
  double x = 0.0 ;
  double y = 0.0 ;

  ierr = 0;
  nn = abs ( n );
/*
  Form L in the arrays B and DL.
*/
  for ( i = 0; i < n; i++ )
  {
    for ( j = i; j < n; j++ )
    {
      x = b[i+j*n];

      for ( k = 0; k < i; k++ )
      {
        x = x - b[i+k*n] * b[j+k*n];
      }

      if ( j == i )
      {
        if ( x <= 0.0 )
        {
          printf ( "\n" );
          printf ( "REDUC - Fatal error!\n" );
          printf ( "   The matrix is not positive definite.\n" );
          ierr = 7 * n + 1;
          return ierr;
        }

        y = sqrt ( x );
        dl[i] = y;
      }
      else
      {
        b[j+i*n] = x / y;
      }
    }
  }
/*
  Form the transpose of the upper triangle of INV(L)*A
  in the lower triangle of the array A.
*/
  for ( i = 0; i < nn; i++ )
  {
    y = dl[i];

    for ( j = i; j < nn; j++ )
    {
      x = a[i+j*n];

      for ( k = 0; k < i; k++ )
      {
        x = x - b[i+k*n] * a[j+k*n];
      }
      a[j+i*n] = x / y;
    }
  }
/*
  Pre-multiply by INV(L) and overwrite.
*/
  for ( j = 0; j < nn; j++ )
  {
    for ( i = j; i < nn; i++ )
    {
      x = a[i+j*n];

      for ( k = j; k < i; k++ )
      {
        x = x - a[k+j*n] * b[i+k*n];
      }

      for ( k = 0; k < j; k++ )
      {
        x = x - a[j+k*n] * b[i+k*n];
      }
      a[i+j*n] = x / dl[i];
    }
  }

  return ierr;
}
/******************************************************************************/

int reduc2 ( int n, double a[], double b[], double dl[] )

/******************************************************************************/
/*
  Purpose:

    REDUC2 reduces the eigenvalue problem A*B*x=lamdba*x to A*x=lambda*x.

  Discussion:

    REDUC2 reduces the generalized symmetric eigenproblems
    a*b*x=lambda*x or b*a*y=lambda*y, where B is positive definite,
    to the standard symmetric eigenproblem using the Cholesky
    factorization of B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrices A and B.  If the
    Cholesky factor L of B is already available, N should be prefixed with a
    minus sign.

    Input/output, double A[N,N].  On input, A contains a real
    symmetric matrix.  Only the full upper triangle of the matrix need be
    supplied.  On output, A contains in its full lower triangle the full lower
    triangle of the symmetric matrix derived from the reduction to the
    standard form.  The strict upper triangle of a is unaltered.

    Input/output, double B[N,N].  On input, the real symmetric
    input matrix.  Only the full upper triangle of the matrix need be supplied.
    If N is negative, the strict lower triangle of B contains, instead, the
    strict lower triangle of its Cholesky factor L.  In any case, on output,
    B contains in its strict lower triangle the strict lower triangle of
    its Cholesky factor L.  The full upper triangle of B is unaltered.

    Input/output, double DL[N].  If N is negative, then the DL
    contains the diagonal elements of L on input.  In any case, DL will contain
    the diagonal elements of L on output,

    Output, int REDUC2, error flag.
    0, for normal return,
    7*N+1, if B is not positive definite.
*/
{
  int i;
  int ierr;
  int j;
  int k;
  int nn;
  double x = 0.0 ;
  double y = 0.0 ;

  ierr = 0;
  nn = abs ( n );
/*
  Form L in the arrays B and DL.
*/
  for ( i = 0; i < n; i++ )
  {
    for ( j = i; j < n; j++ )
    {
      x = b[i+j*nn];

      for ( k = 0; k < i; k++ )
      {
        x = x - b[i+k*nn] * b[j+k*nn];
      }

      if ( j == i )
      {
        if ( x <= 0.0 )
        {
          printf ( "\n" );
          printf ( "REDUC2 - Fatal error!\n" );
          printf ( "  The matrix is not positive definite.\n" );
          ierr = 7 * n + 1;
          return ierr;
        }

        y = sqrt ( x );
        dl[i] = y;
      }
      else
      {
        b[j+i*nn] = x / y;
      }
    }
  }
/*
  Form the lower triangle of A*L in the lower triangle of A.
*/
  for ( i = 0; i < nn; i++ )
  {
    for ( j = 0; j <= i; j++ )
    {
      x = a[j+i*nn] * dl[j];

      for ( k = j + 1; k <= i; k++ )
      {
        x = x + a[k+i*nn] * b[k+j*nn];
      }

      for ( k = i + 1; k < nn; k++ )
      {
        x = x + a[i+k*nn] * b[k+j*nn];
      }
      a[i+j*nn] = x;
    }
  }
/*
  Pre-multiply by L' and overwrite.
*/
  for ( i = 0; i < nn; i++ )
  {
    y = dl[i];

    for ( j = 0; j <= i; j++ )
    {
      x = y * a[i+j*nn];

      for ( k = i + 1; k < nn; k++ )
      {
        x = x + a[k+j*nn] * b[k+i*nn];
      }
      a[i+j*nn] = x;
    }
  }

  return ierr;
}
/******************************************************************************/

int rg_elm ( int n, double a[], double wr[], double wi[], bool matz, 
  double z[] )

/******************************************************************************/
/*
  Purpose:

    RG_ELM computes eigenvalues and eigenvectors of a real general matrix.

  Discussion:

    RG_ELM calls EISPACK routines to find the eigenvalues and eigenvectors
    of a real general matrix using elementary transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double A[N,N], the real general matrix.  On 
    output, A has been overwritten.

    Input, bool MATZ, is false if only eigenvalues are desired, 
    and true if both eigenvalues and eigenvectors are desired.

    Output, double WR[N], WI[N], the real and imaginary parts,
    respectively, of the eigenvalues.  Complex conjugate pairs of eigenvalues
    appear consecutively with the eigenvalue having the positive imaginary
    part first.

    Output, double Z[N,N], contains the real and imaginary parts of 
    the eigenvectors if MATZ is true.  If the J-th eigenvalue is real, the
    J-th column of Z contains its eigenvector.  If the J-th eigenvalue is
    complex with positive imaginary part, the J-th and (J+1)-th columns of
    Z contain the real and imaginary parts of its eigenvector.  The
    conjugate of this vector is the eigenvector for the conjugate eigenvalue.

    Output, int RG_ELM, an error completion code described in 
    the documentation for HQR and HQR2.  The normal completion code is zero.
*/
{
  double *fv1;
  int ierr;
  int is1;
  int is2;
  int *iv1;

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );

  balanc ( n, a, &is1, &is2, fv1 );

  iv1 = ( int * ) malloc ( n * sizeof ( int ) );

  elmhes ( n, is1, is2, a, iv1 );

  if ( ! matz )
  {
    ierr = hqr ( n, is1, is2, a, wr, wi );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( iv1 );
      printf ( "\n" );
      printf ( "RG_ELM - Fatal error!\n" );
      printf ( "  Error return from HQR.\n" );
      return ierr;
    }
  }
  else
  {
    eltran ( n, is1, is2, a, iv1, z );

    ierr = hqr2 ( n, is1, is2, a, wr, wi, z );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( iv1 );
      printf ( "\n" );
      printf ( "RG_ELM - Fatal error!\n" );
      printf ( "  Error return from HQR2.\n" );
      return ierr;
    }

    balbak ( n, is1, is2, fv1, n, z );
  }
/*
  Free memory.
*/
  free ( fv1 );
  free ( iv1 );

  return ierr;
}
/******************************************************************************/

int rg_ort ( int n, double a[], double wr[], double wi[], bool matz, 
  double z[] )

/******************************************************************************/
/*
  Purpose:

    RG_ORT computes eigenvalues and eigenvectors of a real general matrix.

  Discussion:

    RG_ORT calls EISPACK routines to find the eigenvalues and eigenvectors
    of a real general matrix, using orthogonal transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double A[N,N], the real general matrix.  On 
    output, A has been overwritten.

    Input, bool MATZ, is false if only eigenvalues are desired, 
    and true if both eigenvalues and eigenvectors are desired.

    Output, double WR[N], WI[N], the real and imaginary parts,
    respectively, of the eigenvalues.  Complex conjugate pairs of eigenvalues
    appear consecutively with the eigenvalue having the positive imaginary
    part first.

    Output, double Z[N,N], contains the real and imaginary parts of 
    the eigenvectors if MATZ is true.  If the J-th eigenvalue is real, the
    J-th column of Z contains its eigenvector.  If the J-th eigenvalue is
    complex with positive imaginary part, the J-th and (J+1)-th columns of
    Z contain the real and imaginary parts of its eigenvector.  The
    conjugate of this vector is the eigenvector for the conjugate eigenvalue.

    Output, int RG_ORT, an error completion code described in 
    the documentation for HQR and HQR2.  The normal completion code is zero.
*/
{
  double *fv1;
  int ierr;
  int is1;
  int is2;
  double *ort;

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );

  balanc ( n, a, &is1, &is2, fv1 );

  ort = ( double * ) malloc ( n * sizeof ( double ) );

  orthes ( n, is1, is2, a, ort );

  if ( ! matz )
  {
    ierr = hqr ( n, is1, is2, a, wr, wi );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( ort );
      printf ( "\n" );
      printf ( "RG_ORT - Fatal error!\n" );
      printf ( "  Error return from HQR.\n" );
      return ierr;
    }
  }
  else
  {
    ortran ( n, is1, is2, a, ort, z );

    ierr = hqr2 ( n, is1, is2, a, wr, wi, z );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( ort );
      printf ( "\n" );
      printf ( "RG_ORT - Fatal error!\n" );
      printf ( "  Error return from HQR2.\n" );
      return ierr;
    }

    balbak ( n, is1, is2, fv1, n, z );
  }

  free ( fv1 );
  free ( ort );

  return ierr;
}
/******************************************************************************/

int rgg ( int n, double a[], double b[], double alfr[], double alfi[], 
  double beta[], bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    RGG: eigenvalues/vectors for the generalized problem A*x = lambda*B*x.

  Discussion:

    RGG calls the recommended sequence of EISPACK routines
    to find the eigenvalues and eigenvectors (if desired)
    for the real general generalized eigenproblem

      A * x = lambda * B * x.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrices A and B.

    Input/output, double A[N,N], B[N,N], the two real general 
    matrices.  On output, A and B have been overwritten.

    Input, bool MATZ, is false if only eigenvalues are desired, 
    and true if both eigenvalues and eigenvectors are desired.

    Output, double ALFR[N], ALFI[N], the real and imaginary parts,
    respectively, of the numerators of the eigenvalues.

    Output, double BETA[N], the denominators of the eigenvalues,
    which are thus given by the ratios (ALFR + I * ALFI ) / BETA.
    Complex conjugate pairs of eigenvalues appear consecutively
    with the eigenvalue having the positive imaginary part first.

    Output, double Z[N,N], contains the real and imaginary parts of 
    the eigenvectors if MATZ is true.  If the J-th eigenvalue is real, the
    J-th column of Z contains its eigenvector.  If the J-th eigenvalue is
    complex with positive imaginary part, the J-th and (J+1)-th columns of
    Z contain the real and imaginary parts of its eigenvector.  The
    conjugate of this vector is the eigenvector for the conjugate eigenvalue.

    Output, int RGG, is set equal to an error completion code
    described in the documentation for QZIT.  The normal completion
    code is zero.
*/
{
  double eps1;
  int ierr;

  eps1 = 0.0;

  qzhes ( n, a, b, matz, z );

  ierr = qzit ( n, a, b, eps1, matz, z );

  if ( ierr != 0 )
  {
    printf ( "\n" );
    printf ( "RGG - Fatal error!\n" );
    printf ( "  Error return from QZIT.\n" );
    return ierr;
  }

  qzval ( n, a, b, alfr, alfi, beta, matz, z );

  if ( matz ) 
  {
    qzvec ( n, a, b, alfr, alfi, beta, z );
  }

  return ierr;
}
/******************************************************************************/

int rs ( int n, double a[], double w[], bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    RS computes eigenvalues and eigenvectors of real symmetric matrix.

  Discussion:

    RS calls the recommended sequence of
    function from the eigensystem subroutine package (eispack)
    to find the eigenvalues and eigenvectors (if desired)
    of a real symmetric matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the real symmetric matrix.

    Input, bool MATZ, is false if only eigenvalues are desired, 
    and true if both eigenvalues and eigenvectors are desired.

    Output, double W[N], the eigenvalues in ascending order.

    Output, double Z[N*N], contains the eigenvectors, if MATZ
    is true.

    Output, int RS, is set equal to an error
    completion code described in the documentation for TQLRAT and TQL2.
    The normal completion code is zero.
*/
{
  double *fv1;
  double *fv2;
  int ierr;

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );

  if ( ! matz )
  {
    fv2 = ( double * ) malloc ( n * sizeof ( double ) );

    tred1 ( n, a, w, fv1, fv2 );

    ierr = tqlrat ( n, w, fv2 );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );

      printf ( "\n" );
      printf ( "RS - Fatal error!\n" );
      printf ( "  Error return from TQLRAT!\n" );
      return ierr;
    }

    free ( fv2 );
  }
  else
  {
    tred2 ( n, a, w, fv1, z );

    ierr = tql2 ( n, w, fv1, z );

    if ( ierr != 0 )
    {
      free ( fv1 );

      printf ( "\n" );
      printf ( "RS - Fatal error!\n" );
      printf ( "  Error return from TQL2!\n" );
      return ierr;
    }
  }

  free ( fv1 );

  return ierr;
}
/******************************************************************************/

int rsb ( int n, int mb, double a[], double w[], bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    RSB computes eigenvalues and eigenvectors of a real symmetric band matrix.

  Discussion:

    RSB calls the recommended sequence of
    functions from the eigensystem subroutine package (eispack)
    to find the eigenvalues and eigenvectors (if desired)
    of a real symmetric band matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int MB, the half band width of the matrix,
    defined as the number of adjacent diagonals, including the principal
    diagonal, required to specify the non-zero portion of the lower triangle
    of the matrix.

    Input, double A[N*MB], contains the lower triangle of the real
    symmetric band matrix.  Its lowest subdiagonal is stored in the last N+1-MB
    positions of the first column, its next subdiagonal in the last
    N+2-MB positions of the second column, further subdiagonals similarly,
    and finally its principal diagonal in the N positions of the last
    column.  Contents of storages not part of the matrix are arbitrary.

    Input, bool MATZ, is false if only eigenvalues are desired,
    and true if both eigenvalues and eigenvectors are desired.

    Output, double W[N], the eigenvalues in ascending order.

    Output, double Z[N*N], contains the eigenvectors, if MATZ
    is true.

    Output, int BANDR, is set to an error
    completion code described in the documentation for TQLRAT and TQL2.
    The normal completion code is zero.
*/
{
  double *fv1;
  double *fv2;
  int ierr;

  if ( mb <= 0 )
  {
    printf ( "\n" );
    printf ( "RSB - Fatal error!\n" );
    printf ( "  MB <= 0!\n" );
    ierr = 12 * n;
    return ierr;
  }

  if ( n < mb )
  {
    printf ( "\n" );
    printf ( "RSB - Fatal error!\n" );
    printf ( "  N < MB!\n" );
    ierr = 12 * n;
    return ierr;
  }

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );

  if ( ! matz )
  {
    fv2 = ( double * ) malloc ( n * sizeof ( double ) );

    bandr ( n, mb, a, w, fv1, fv2, matz, z );

    ierr = tqlrat ( n, w, fv2 );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );
      printf ( "\n" );
      printf ( "RSB - Fatal error!\n" );
      printf ( "  Error return from TQLRAT!\n" );
      return ierr;
    }

    free ( fv2 );
  }
  else
  {
    bandr ( n, mb, a, w, fv1, fv1, matz, z );

    ierr = tql2 ( n, w, fv1, z );

    if ( ierr != 0 )
    {
      free ( fv1 );

      printf ( "\n" );
      printf ( "RSB - Fatal error!\n" );
      printf ( "  Error return from TQLRAT!\n" );
      return ierr;
    }
  }

  free ( fv1 );

  return ierr;
}
/******************************************************************************/

int rsg ( int n, double a[], double b[], double w[], bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    RSG computes eigenvalues/vectors, A*x=lambda*B*x, A symmetric, B pos-def.

  Discussion:

    RSG calls the recommended sequence of EISPACK functions
    to find the eigenvalues and eigenvectors (if desired)
    for the real symmetric generalized eigenproblem  a x = (lambda) b x.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrices A and B.

    Input, double A[N,N], contains a real symmetric matrix.

    Input, double B[N,N], contains a positive definite real
    symmetric matrix.

    Input, bool MATZ, is false if only eigenvalues are desired, 
    and true if both eigenvalues and eigenvectors are desired.

    Output, double W[N], the eigenvalues in ascending order.

    Output, double Z[N,N], contains the eigenvectors, if MATZ
    is true.

    Output, int RSG, is set to an error
    completion code described in the documentation for TQLRAT and TQL2.
    The normal completion code is zero.
*/
{
  double *fv1;
  double *fv2;
  int ierr;

  fv2 = ( double * ) malloc ( n * sizeof ( double ) );

  ierr = reduc ( n, a, b, fv2 );

  if ( ierr != 0 )
  {
    free ( fv2 );
    printf ( "\n" );
    printf ( "RSG - Fatal error!\n" );
    printf ( "  Error return from REDUC.\n" );
    return ierr;
  }

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );

  if ( ! matz )
  {
    tred1 ( n, a, w, fv1, fv2 );

    ierr = tqlrat ( n, w, fv2 );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );
      printf ( "\n" );
      printf ( "RSG - Fatal error!\n" );
      printf ( "  Error return from TQLRAT!\n" );
      return ierr;
    }
  }
  else
  {
    tred2 ( n, a, w, fv1, z );

    ierr = tql2 ( n, w, fv1, z );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );
      printf ( "\n" );
      printf ( "RSG - Fatal error!\n" );
      printf ( "  Error return from TQL2!\n" );
      return ierr;
    }

    rebak ( n, b, fv2, n, z );
  }

  free ( fv1 );
  free ( fv2 );

  return ierr;
}
/******************************************************************************/

int rsgab ( int n, double a[], double b[], double w[], bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    RSGAB computes eigenvalues/vectors, A*B*x=lambda*x, A symmetric, B pos-def.

  Discussion:

    RSGAB calls the recommended sequence of EISPACK functions
    to find the eigenvalues and eigenvectors (if desired)
    for the real symmetric generalized eigenproblem  
      A * B * x = lambda * x.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrices A and B.

    Input, double A[N,N], contains a real symmetric matrix.

    Input, double B[N,N], contains a positive definite real
    symmetric matrix.

    Input, bool MATZ, is false if only eigenvalues are desired, 
    and true if both eigenvalues and eigenvectors are desired.

    Output, double W[N], the eigenvalues in ascending order.

    Output, double Z[N,N], contains the eigenvectors, if MATZ 
    is true.

    Output, int RSGAB, is set to an error
    completion code described in the documentation for TQLRAT and TQL2.
    The normal completion code is zero.
*/
{
  double *fv1;
  double *fv2;
  int ierr;

  fv2 = ( double * ) malloc ( n * sizeof ( double ) );

  ierr = reduc2 ( n, a, b, fv2 );

  if ( ierr != 0 )
  {
    free ( fv2 );
    printf ( "\n" );
    printf ( "RSGAB - Fatal error!\n" );
    printf ( "  Error return from REDUC2!\n" );
    return ierr;
  }

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );

  if ( ! matz )
  {
    tred1 ( n, a, w, fv1, fv2 );

    ierr = tqlrat ( n, w, fv2 );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );
      printf ( "\n" );
      printf ( "RSGAB - Fatal error!\n" );
      printf ( "  Error return from TQLRAT!\n" );
      return ierr;
    }
  }
  else
  {
    tred2 ( n, a, w, fv1, z );

    ierr = tql2 ( n, w, fv1, z );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );
      printf ( "\n" );
      printf ( "RSB - Fatal error!\n" );
      printf ( "  Error return from TQL2!\n" );
      return ierr;
    }

    rebak ( n, b, fv2, n, z );
  }
/*
  Free memory.
*/
  free ( fv1 );
  free ( fv2 );

  return ierr;
}
/******************************************************************************/

int rsgba ( int n, double a[], double b[], double w[], bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    RSGBA computes eigenvalues/vectors, B*A*x=lambda*x, A symmetric, B pos-def.

  Discussion:

    RSGBA calls the recommended sequence of EISPACK functions
    to find the eigenvalues and eigenvectors (if desired)
    for the real symmetric generalized eigenproblem:

      B * A * x = lambda * x

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrices A and B.

    Input, double A[N,N], a real symmetric matrix.

    Input, double B[N,N], a positive definite symmetric matrix.

    Input, bool MATZ, is false if only eigenvalues are desired, 
    and true if both eigenvalues and eigenvectors are desired.

    Output, double W[N], the eigenvalues in ascending order.

    Output, double Z[N,N], contains the eigenvectors, if MATZ
    is true.

    Output, int RSGBA, is set to an error
    completion code described in the documentation for TQLRAT and TQL2.
    The normal completion code is zero.
*/
{
  double *fv1;
  double *fv2;
  int ierr;

  fv2 = ( double * ) malloc ( n * sizeof ( double ) );

  ierr = reduc2 ( n, a, b, fv2 );

  if ( ierr != 0 )
  {
    free ( fv2 );
    printf ( "\n" );
    printf ( "RSGBA - Fatal error!\n" );
    printf ( "  Error return from REDUC2!\n" );
    return ierr;
  }

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );

  if ( ! matz )
  {
    tred1 ( n, a, w, fv1, fv2 );

    ierr = tqlrat ( n, w, fv2 );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );
      printf ( "\n" );
      printf ( "RSGBA - Fatal error!\n" );
      printf ( "  Error return from TQLRAT!\n" );
      return ierr;
    }
  }
  else
  {
    tred2 ( n, a, w, fv1, z );

    ierr = tql2 ( n, w, fv1, z );

    if ( ierr != 0 )
    {
      free ( fv1 );
      free ( fv2 );
      printf ( "\n" );
      printf ( "RSGBA - Fatal error!\n" );
      printf ( "  Error return from TQL2!\n" );
      return ierr;
    }

    rebakb ( n, b, fv2, n, z );
  }
/*
  Free memory.
*/
  free ( fv1 );
  free ( fv2 );

  return ierr;
}
/******************************************************************************/

int rsm ( int n, double a[], double w[], int m, double z[] )

/******************************************************************************/
/*
  Purpose:

    RSM computes eigenvalues, some eigenvectors, real symmetric matrix.

  Discussion:

    RSM calls the recommended sequence of EISPACK routines
    to find all of the eigenvalues and some of the eigenvectors
    of a real symmetric matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N,N], the symmetric matrix.

    Output, double W[N], the eigenvalues in ascending order.

    Input, int M, specifies the number of eigenvectors to 
    compute.

    Output, double Z[N,M], contains the orthonormal eigenvectors
    associated with the first M eigenvalues.

    Output, int RSM, is set to an error
    completion code described in the documentation for TQLRAT, IMTQLV and
    TINVIT.  The normal completion code is zero.
*/
{
  double *fwork1;
  double *fwork2;
  double *fwork3;
  int ierr;
  int *iwork;

  fwork1 = ( double * ) malloc ( n * sizeof ( double ) );
  fwork2 = ( double * ) malloc ( n * sizeof ( double ) );

  if ( m <= 0 )
  {
    tred1 ( n, a, w, fwork1, fwork2 );

    ierr = tqlrat ( n, w, fwork2 );

    if ( ierr != 0 )
    {
      free ( fwork1 );
      free ( fwork2 );
      printf ( "\n" );
      printf ( "RSM - Fatal error!\n" );
      printf ( "  Error return from TQLRAT!\n" );
      return ierr;
    }
  }
  else
  {
    fwork3 = ( double * ) malloc ( n * sizeof ( double ) );

    tred1 ( n, a, fwork1, fwork2, fwork3 );

    iwork = ( int * ) malloc ( n * sizeof ( int ) );

    ierr = imtqlv ( n, fwork1, fwork2, fwork3, w, iwork );

    if ( ierr != 0 )
    {
      free ( fwork1 );
      free ( fwork2 );
      free ( fwork3 );
      free ( iwork );
      printf ( "\n" );
      printf ( "RSM - Fatal error!\n" );
      printf ( "  Error return from IMTQLV!\n" );
      return ierr;
    }

    ierr = tinvit ( n, fwork1, fwork2, fwork3, m, w, iwork, z );

    if ( ierr != 0 )
    {
      free ( fwork1 );
      free ( fwork2 );
      free ( fwork3 );
      free ( iwork );
      printf ( "\n" );
      printf ( "RSM - Fatal error!\n" );
      printf ( "  Error return from TINVIT!\n" );
      return ierr;
    }

    trbak1 ( n, a, fwork2, m, z );

    free ( fwork3 );
    free ( iwork );
  }

  free ( fwork1 );
  free ( fwork2 );

  return ierr;
}
/******************************************************************************/

int rsp ( int n, int nv, double a[], double w[], bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    RSP computes eigenvalues and eigenvectors of real symmetric packed matrix.

  Discussion:

    RSP calls the recommended sequence of EISPACK routines
    to find the eigenvalues and eigenvectors (if desired)
    of a real symmetric packed matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int NV, the dimension of the array A, which
    must be at least (N*(N+1))/2.

    Input, double A[NV], contains the lower triangle of the
    real symmetric packed matrix stored row-wise.

    Input, bool MATZ, is false if only eigenvalues are desired,
    and true if both eigenvalues and eigenvectors are desired.

    Output, double W[N], the eigenvalues in ascending order.

    Output, double Z[N,N], contains the eigenvectors, if MATZ is 
    true.

    Output, int RSP, is set to an error
    completion code described in the documentation for TQLRAT and TQL2.
    The normal completion code is zero.
*/
{
  double *fv1;
  double *fv2;
  int ierr;

  if ( ( n * ( n + 1 ) ) / 2 > nv )
  {
    ierr = 20 * n;
    printf ( "\n" );
    printf ( "RSP - Fatal error!\n" );
    printf ( "  NV is too small!\n" );
    return ierr;
  }

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );
  fv2 = ( double * ) malloc ( n * sizeof ( double ) );

  tred3 ( n, nv, a, w, fv1, fv2 );

  if ( ! matz )
  {
    ierr = tqlrat ( n, w, fv2 );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "RSP - Fatal error!\n" );
      printf ( "  Error return from TQLRAT.\n" );
      free ( fv1 );
      free ( fv2 );
      return ierr;
    }
  }
  else
  {
    r8mat_identity ( n, z );

    ierr = tql2 ( n, w, fv1, z );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "RSP - Fatal error!\n" );
      printf ( "  Error return from TQL2.\n" );
      free ( fv1 );
      free ( fv2 );
      return ierr;
    }

    trbak3 ( n, nv, a, n, z );
  }

  free ( fv1 );
  free ( fv2 );

  return ierr;
}
/******************************************************************************/

int rspp ( int n, int nv, double a[], double w[], bool matz, double z[], int m, 
  bool type )

/******************************************************************************/
/*
  Purpose:

    RSPP computes some eigenvalues/vectors, real symmetric packed matrix.

  Discussion:

    RSPP calls the appropriate routines for the following problem:

    Given a symmetric matrix A, which is stored in a packed mode, find
    the M smallest or largest eigenvalues, and corresponding eigenvectors.

    The routine RSP returns all eigenvalues and eigenvectors.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of A, the number of rows and
    columns in the original matrix.

    Input, int NV, is the of the array A as specified in the
    calling program.  NV must not be less than N*(N+1)/2.

    Input, double A[(N*(N+1))/2], on input the lower triangle of the
    real symmetric matrix, stored row-wise in the vector,
    in the order A(1,1), / A(2,1), A(2,2), / A(3,1), A(3,2), A(3,3)/
    and so on.

    Output, double W[M], the eigenvalues requested.

    Input, bool MATZ, is false if only eigenvalues are desired,
    and true if both eigenvalues and eigenvectors are desired.

    Output, double Z[N,M], the eigenvectors.

    Input, int M, the number of eigenvalues to be found.

    Input, bool TYPE, is true if the smallest eigenvalues
    are to be found, or false if the largest ones are sought.

    Output, int RSPP, error flag from RATQR.  
    0 on normal return. 
    nonzero, the algorithm broke down while computing an eigenvalue.
*/
{
  double *bd;
  double eps1;
  int idef;
  int ierr;
  int *iwork;
  double *work1;
  double *work2;
  double *work3;
/*
  IDEF =
    -1 if the matrix is known to be negative definite,
    +1 if the matrix is known to be positive definite, or
    0 otherwise.
*/
  idef = 0;
/*
  Reduce to symmetric tridiagonal form.
*/
  work1 = ( double * ) malloc ( n * sizeof ( double ) );
  work2 = ( double * ) malloc ( n * sizeof ( double ) );
  work3 = ( double * ) malloc ( n * sizeof ( double ) );

  tred3 ( n, nv, a, work1, work2, work3 );
/*
  Find the eigenvalues.
*/
  eps1 = 0.0;

  iwork = ( int * ) malloc ( n * sizeof ( int ) );
  bd = ( double * ) malloc ( n * sizeof ( double ) );

  ierr = ratqr ( n, eps1, work1, work2, work3, m, w, iwork, bd, type, idef );

  if ( ierr != 0 )
  {
    free ( bd );
    free ( iwork );
    free ( work1 );
    free ( work2 );
    free ( work3 );
    printf ( "\n" );
    printf ( "RSPP - Fatal error!\n" );
    printf ( "  Error return from RATQR.\n" );
    return ierr;
  }
/*
  Find eigenvectors for the first M eigenvalues.
*/
  if ( matz )
  {
    ierr = tinvit ( n, work1, work2, work3, m, w, iwork, z );

    if ( ierr != 0 )
    {
      free ( bd );
      free ( iwork );
      free ( work1 );
      free ( work2 );
      free ( work3 );
      printf ( "\n" );
      printf ( "RSPP - Fatal error!\n" );
      printf ( "  Error return from TINVIT.\n" );
      return ierr;
    }
/*
  Reverse the transformation.
*/
    trbak3 ( n, nv, a, m, z );
  }

  free ( bd );
  free ( iwork );
  free ( work1 );
  free ( work2 );
  free ( work3 );

  return ierr;
}
/******************************************************************************/

int rst ( int n, double w[], double e[], bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    RST computes eigenvalues/vectors, real symmetric tridiagonal matrix.

  Discussion:

    RST calls the recommended sequence of EISPACK routines
    to find the eigenvalues and eigenvectors (if desired)
    of a real symmetric tridiagonal matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double W[N].  On input, the diagonal elements
    of the real symmetric tridiagonal matrix.  On output, the eigenvalues in
    ascending order.

    Input, double E[N], the subdiagonal elements of the matrix in
    E(2:N).  E(1) is arbitrary.

    Input, bool MATZ, is false if only eigenvalues are desired,
    and true if both eigenvalues and eigenvectors are desired.

    Output, double Z[N,N], contains the eigenvectors, if MATZ
    is true.

    Output, int RST_TEST, is set to an error
    completion code described in the documentation for IMTQL1 and IMTQL2.
    The normal completion code is zero.
*/
{
  int ierr;

  if ( ! matz )
  {
    ierr = imtql1 ( n, w, e );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "RST - Fatal error!\n" );
      printf ( "  Error return from IMTQL1.\n" );
      return ierr;
    }
  }
  else
  {
    r8mat_identity ( n, z );

    ierr = imtql2 ( n, w, e, z );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "RST - Fatal error!\n" );
      printf ( "  Error return from IMTQL2.\n" );
      return ierr;
    }
  }

  return ierr;
}
/******************************************************************************/

int rt ( int n, double a[], double w[], bool matz, double z[] )

/******************************************************************************/
/*
  Purpose:

    RT computes eigenvalues/vectors, real sign-symmetric tridiagonal matrix.

  Discussion:

    RT calls the recommended sequence of EISPACK routines
    to find the eigenvalues and eigenvectors (if desired)
    of a special real tridiagonal matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N,N], contains the special real tridiagonal
    matrix in its first three columns.  The subdiagonal elements are stored
    in the last N-1 positions of the first column, the diagonal elements
    in the second column, and the superdiagonal elements in the first N-1
    positions of the third column.  Elements A(1,1) and A(N,3) are arbitrary.

    Input, bool MATZ, is false if only eigenvalues are desired,
    and true if both eigenvalues and eigenvectors are desired.

    Output, double W[N], the eigenvalues in ascending order.

    Output, double Z[N,N], contains the eigenvectors, if MATZ is true.

    Output, int RT, is set to an error
    completion code described in the documentation for IMTQL1 and IMTQL2.
    The normal completion code is zero.
*/
{
  double *fv1;
  int ierr;

  fv1 = ( double * ) malloc ( n * sizeof ( double ) );

  if ( ! matz )
  {
    ierr = figi ( n, a, w, fv1, fv1 );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "RT - Fatal error!\n" );
      printf ( "  Error return from FIGI.\n" );
      return ierr;
    }

    ierr = imtql1 ( n, w, fv1 );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "RT - Fatal error!\n" );
      printf ( "  Error return from IMTQL1.\n" );
      return ierr;
    }
  }
  else
  {
    ierr = figi2 ( n, a, w, fv1, z );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "RT - Fatal error!\n" );
      printf ( "  Error return from FIGI2.\n" );
      return ierr;
    }

    ierr = imtql2 ( n, w, fv1, z );

    if ( ierr != 0 )
    {
      printf ( "\n" );
      printf ( "RT - Fatal error!\n" );
      printf ( "  Error return from IMTQL2.\n" );
      return ierr;
    }
  }

  free ( fv1 );

  return ierr;
}
/******************************************************************************/

int sturm_sequence ( double d[], double e[], double e2[], int n, int p, int q, 
  double x1 )

/******************************************************************************/
/*
  Purpose:

    STURM_SEQUENCE counts eigenvalues of a symmetric tridiagonal submatrix.

  Discussion:

    Let A be a symmetric tridiagonal matrix, and consider the submatrix
    defined by indices P through Q.

    STURM_SEQUENCE will determine the number of eigenvalues associated
    with this submatrix which are no larger than a given upper bound X1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, double D[N], the diagonal elements of the input matrix.

    Input, double E[N], the subdiagonal elements of the matrix.  
    E(1) is arbitrary.

    Input, double E2[N], the squares of the corresponding
    elements of E, with zeros corresponding to negligible elements of E.

    Input, int N, the order of the matrix.

    Input, int P, Q, the lower and upper limits on the submatrix.
    Note that these are 0-based indices.  0 <= P <= Q < N.

    Input, double X1, an upper bound.

    Output, int STURM_SEQUENCE, the number of submatrix eigenvalues
    less than or equal to X1.
*/
{
  int i;
  int s;
  double u;
  double v;

  s = p;
  u = 1.0;

  for ( i = p; i <= q; i++ )
  {
    if ( u == 0.0 )
    {
      v = fabs ( e[i] ) / r8_epsilon ( );
      if ( e2[i] == 0.0 )
      {
        v = 0.0;
      }
    }
    else
    {
      v = e2[i] / u;
    }

    u = d[i] - x1 - v;
    if ( u < 0.0 )
    {
      s = s + 1;
    }
  }

  return s;
}
/******************************************************************************/

int svd ( int m, int n, double a[], double w[], bool matu, double u[], 
  bool matv, double v[] )

/******************************************************************************/
/*
  Purpose:

    SVD computes the singular value decomposition for a real matrix.

  Discussion:

    SVD determines the singular value decomposition

      A = U * S * V'

    of a real M by N rectangular matrix.  Householder bidiagonalization
    and a variant of the QR algorithm are used.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Golub, Christian Reinsch,
    Numerische Mathematik,
    Volume 14, 1970, pages 403-420.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int M, the number of rows of A and U.

    Input, int N, the number of columns of A and U, and
    the order of V.

    Input, double A[M,N], the M by N matrix to be decomposed.

    Output, double W[N], the singular values of A.  These are the
    diagonal elements of S.  They are unordered.  If an error break; is
    made, the singular values should be correct for indices
    IERR+1, IERR+2,..., N.

    Input, bool MATU, should be set to TRUE if the U matrix in the
    decomposition is desired, and to FALSE otherwise.

    Output, double U[M,N], contains the matrix U, with orthogonal
    columns, of the decomposition, if MATU has been set to TRUE.  Otherwise
    U is used as a temporary array.  U may coincide with A.
    If an error break; is made, the columns of U corresponding
    to indices of correct singular values should be correct.

    Input, bool MATV, should be set to TRUE if the V matrix in the
    decomposition is desired, and to FALSE otherwise.

    Output, double V[N,N], the orthogonal matrix V of the
    decomposition if MATV has been set to TRUE.  Otherwise V is not referenced.
    V may also coincide with A if U is not needed.  If an error
    break; is made, the columns of V corresponding to indices of
    correct singular values should be correct.

    Output, int SVD, error flag.
    0, for normal return,
    K, if the K-th singular value has not been determined after 30 iterations.
*/
{
  double c;
  double f;
  double g;
  double h;
  int i;
  int ierr;
  int ii;
  int its;
  int j;
  int k;
  int l;
  int mn;
  double *rv1;
  double s;
  double scale;
  bool skip;
  double tst1;
  double tst2;
  double x;
  double y;
  double z;

  rv1 = ( double * ) malloc ( n * sizeof ( double ) );

  ierr = 0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      u[i+j*m] = a[i+j*m];
    }
  }
/*
  Householder reduction to bidiagonal form.
*/
  g = 0.0;
  scale = 0.0;
  x = 0.0;

  for ( i = 0; i < n; i++ )
  {
    rv1[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;

    if ( i <= m )
    {
      scale = 0.0;
      for ( ii = i; ii < m; ii++ )
      {
        scale = scale + fabs ( u[ii+i*m] );
      }

      if ( scale != 0.0 )
      {
        for ( ii = i; ii < m; ii++ )
        {
          u[ii+i*m] = u[ii+i*m] / scale;
        }
        s = 0.0;
        for ( ii = i; ii < m; ii++ )
        {
          s = s + u[ii+i*m] * u[ii+i*m];
        }
        f = u[i+i*m];
        g = - sqrt ( s ) * r8_sign ( f );
        h = f * g - s;
        u[i+i*m] = f - g;

        for ( j = i + 1; j < n; j++ )
        {
          s = 0.0;
          for ( k = i; k < m; k++ )
          {
            s = s + u[k+i*m] * u[k+j*m];
          }
          for ( k = i; k < m; k++ )
          {
            u[k+j*m] = u[k+j*m] + s * u[k+i*m] / h;
          }
        }
        for ( k = i; k < m; k++ )
        {
          u[k+i*m] = scale * u[k+i*m];
        }
      }
    }

    w[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;

    if ( i < m && i != n - 1 )
    {
      scale = 0.0;
      for ( k = i + 1; k < n; k++ )
      {
        scale = scale + fabs ( u[i+k*m] );
      }
      if ( scale != 0.0 )
      {
        for ( k = i + 1; k < n; k++ )
        {
          u[i+k*m] = u[i+k*m] / scale;
        }
        s = 0.0;
        for ( k = i + 1; k < n; k++ )
        {
          s = s + u[i+k*m] * u[i+k*m];
        }
        f = u[i+(i+1)*m];
        g = - sqrt ( s ) * r8_sign ( f );
        h = f * g - s;
        u[i+(i+1)*m] = f - g;
        for ( ii = i + 1; ii < n; ii++ )
        {
          rv1[ii] = u[i+ii*m] / h;
        }
        for ( j = i + 1; j < m; j++ )
        {
          s = 0.0;
          for ( k = i + 1; k < n; k++ )
          {
            s = s + u[j+k*m] * u[i+k*m];
          }
          for ( k = i + 1; k < n; k++ )
          {
            u[j+k*m] = u[j+k*m] + s * rv1[k];
          }
        }
        for ( k = i + 1; k < n; k++ )
        {
          u[i+k*m] = u[i+k*m] * scale;
        }
      }
    }

    x = r8_max ( x, fabs ( w[i] ) + fabs ( rv1[i] ) );
  }
/*
  Accumulation of right-hand transformations.
*/
  if ( matv )
  {
    for ( i = n - 1; 0 <= i; i-- )
    {
      if ( i < n - 1 )
      {
        if ( g != 0.0 )
        {
          for ( k = i + 1; k < n; k++ )
          {
            v[k+i*n] = ( u[i+k*m] / u[i+(i+1)*m] ) / g;
          }
          for ( j = i + 1; j < n; j++ )
          {
            s = 0.0;
            for ( k = i + 1; k < n; k++ )
            {
              s = s + u[i+k*m] * v[k+j*n];
            }
            for ( k = i + 1; k < n; k++ )
            {
              v[k+j*n] = v[k+j*n] + s * v[k+i*n];
            }
          }
        }
        for ( k = i + 1; k < n; k++ )
        {
          v[i+k*n] = 0.0;
          v[k+i*n] = 0.0;
        }
      }
      v[i+i*n] = 1.0;
      g = rv1[i];
    }
  }
/*
  Accumulation of left-hand transformations.
*/
  if ( matu )
  {
    mn = i4_min ( m, n );

    for ( i = mn - 1; 0 <= i; i-- )
    {
      g = w[i];

      if ( i != n - 1 )
      {
        for ( k = i + 1; k < n; k++ )
        {
          u[i+k*m] = 0.0;
        }
      }

      if ( g != 0.0 )
      {
        if ( i != mn - 1 )
        {
          for ( j = i + 1; j < n; j++ )
          {
            s = 0.0;
            for ( k = i + 1; k < m; k++ )
            {
              s = s + u[k+i*m] * u[k+j*m];
            }
            f = ( s / u[i+i*m] ) / g;
            for ( k = i; k < m; k++ )
            {
              u[k+j*m] = u[k+j*m] + f * u[k+i*m];
            }
          }
        }
        for ( k = i; k < m; k++ )
        {
          u[k+i*m] = u[k+i*m] / g;
        }
      }
      else
      {
        for ( k = i; k < m; k++ )
        {
          u[k+i*m] = 0.0;
        }
      }
      u[i+i*m] = u[i+i*m] + 1.0;
    }
  }
/*
  Diagonalization of the bidiagonal form.
*/
  tst1 = x;

  for ( k = n - 1; 0 <= k; k-- )
  {
    its = 0;
/*
  Test for splitting.
*/
    while ( true )
    {
      skip = false;

      for ( l = k; 1 <= l; l-- )
      {
        tst2 = tst1 + fabs ( rv1[l] );

        if ( tst2 == tst1 )
        {
          skip = true;
          break;
        }

        tst2 = tst1 + fabs ( w[l-1] );

        if ( tst2 == tst1 )
        {
          break;
        }
      }
/*
  Cancellation of rv1[l] if L greater than 1.
*/
      if ( ! skip )
      {
        c = 0.0;
        s = 1.0;

        for ( i = l; i <= k; i++ )
        {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];
          tst2 = tst1 + fabs ( f );

          if ( tst2 == tst1 )
          {
            break;
          }

          g = w[i];
          h = pythag ( f, g );
          w[i] = h;
          c = g / h;
          s = - f / h;

          if ( matu )
          {
            for ( j = 0; j < m; j++ )
            {
              y = u[j+(l-1)*m];
              z = u[j+i*m];
              u[j+(l-1)*m] =  y * c + z * s;
              u[j+i*m] =    - y * s + z * c;
            }
          }
        }
      }

      z = w[k];
/*
  Convergence.
*/
      if ( l == k )
      {
        if ( z <= 0.0 )
        {
          w[k] = - z;
          if ( matv )
          {
            for ( j = 0; j < n; j++ )
            {
              v[j+k*n] = - v[j+k*n];
            }
          }
        }
        break;
      }
/*
  Shift from bottom 2 by 2 minor.
*/
      if ( 30 <= its )
      {
        free ( rv1 );
        ierr = k;
        return ierr;
      }

      its = its + 1;
      x = w[l];
      y = w[k-1];
      g = rv1[k-1];
      h = rv1[k];
      f = 0.5 * ( ( ( g + z ) / h ) * ( ( g - z ) / y ) + y / h - h / y );
      g = pythag ( f, 1.0 );
      f = x - ( z / x ) * z + ( h / x ) 
        * ( y / ( f + fabs ( g ) * r8_sign ( f ) ) - h );
/*
  Next QR transformation.
*/
      c = 1.0;
      s = 1.0;

      for ( i = l + 1; i <= k; i++ )
      {
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = pythag ( f, h );
        rv1[i-1] = z;
        c = f / z;
        s = h / z;
        f =   x * c + g * s;
        g = - x * s + g * c;
        h = y * s;
        y = y * c;

        if ( matv )
        {
          for ( j = 0; j < n; j++ )
          {
            x = v[j+(i-1)*n];
            z = v[j+i*n];
            v[j+(i-1)*n] =   x * c + z * s;
            v[j+i*n]     = - x * s + z * c;
          }
        }

        z = pythag ( f, h );
        w[i-1] = z;
/*
  Rotation can be arbitrary if Z is zero.
*/
        if ( z != 0.0 )
        {
          c = f / z;
          s = h / z;
        }

        f =   c * g + s * y;
        x = - s * g + c * y;

        if ( matu )
        {
          for ( j = 0; j < m; j++ )
          {
            y = u[j+(i-1)*m];
            z = u[j+i*m];
            u[j+(i-1)*m] =   y * c + z * s;
            u[j+i*m]     = - y * s + z * c;
          }
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }

  free ( rv1 );

  return ierr;
}
/******************************************************************************/

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

int tinvit ( int n, double d[], double e[], double e2[], int m, double w[], 
  int ind[], double z[] )

/******************************************************************************/
/*
  Purpose:

    TINVIT computes eigenvectors from eigenvalues, real tridiagonal symmetric.

  Discussion:

    TINVIT finds eigenvectors of a tridiagonal symmetric matrix corresponding 
    to specified eigenvalues using inverse iteration.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 February 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, 
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double D[N], the diagonal elements of the matrix.

    Input, double E[N], contains the subdiagonal elements of
    the input matrix in E(2:N).  E(1) is arbitrary.

    Input, double E2[N], contains the squares of the corresponding
    elements of E, with zeros corresponding to negligible elements of E.
    E(I) is considered negligible if it is not larger than the product of
    the relative machine precision and the sum of the magnitudes of D(I)
    and D(I-1).  E2(1) must contain 0.0 if the eigenvalues are in
    ascending order, or 2.0 if the eigenvalues are in descending order.
    If BISECT, TRIDIB, or IMTQLV has been used to find the eigenvalues,
    their output E2 array is exactly what is expected here.

    Input, int M, the number of specified eigenvalues.

    Input, double W[M], the eigenvalues.

    Input, int IND[M], the submatrix indices associated with 
    the corresponding eigenvalues in W: 1 for eigenvalues belonging to the
    first submatrix from the top, 2 for those belonging to the second
    submatrix, and so on.

    Output, double Z[N,M], the associated set of orthonormal
    eigenvectors.  Any vector which fails to converge is set to zero.

    Output, int TINVIT, error flag.
    0, for normal return,
    -R, if the eigenvector corresponding to the R-th eigenvalue fails to
      converge in 5 iterations.
*/
{ 
  double eps2 = 0.0 ;
  double eps3 = 0.0 ;
  double eps4 = 0.0 ;
  int group = 0;
  int i;
  int ierr;
  int its;
  int j;
  int jj;
  int k;
  double norm;
  double order;
  int p;
  int q;
  int r;
  double *rv1;
  double *rv2;
  double *rv3;
  double *rv4;
  double *rv6;
  int s;
  int tag;
  double u;
  double uk = 0.0 ;
  double v;
  double x0;
  double x1;
  double xu = 0.0 ;

  ierr = 0;

  if ( m == 0 )
  {
    return ierr;
  }

  rv1 = ( double * ) malloc ( n * sizeof ( double ) );
  rv2 = ( double * ) malloc ( n * sizeof ( double ) );
  rv3 = ( double * ) malloc ( n * sizeof ( double ) );
  rv4 = ( double * ) malloc ( n * sizeof ( double ) );
  rv6 = ( double * ) malloc ( n * sizeof ( double ) );

  u = 0.0;
  x0 = 0.0;

  tag = -1;
  order = 1.0 - e2[0];
  q = -1;
/*
  Establish and process next submatrix.
*/
  while ( true )
  {
    p = q + 1;

    for ( q = p; q < n; q++ )
    {
      if ( q == n - 1 )
      {
        break;
      }
      if ( e2[q+1] == 0.0 )
      {
        break;
      }
    }
/*
  Find vectors by inverse iteration.
*/
    tag = tag + 1;
    s = 0;

    for ( r = 0; r < m; r++ )
    {
      if ( ind[r] != tag )
      {
        continue;
      }

      its = 1;
      x1 = w[r];
/*
  Look for close or coincident roots.
*/
      if ( s != 0 )
      {
        if ( eps2 <= fabs ( x1 - x0 ) )
        {
          group = 0;
        }
        else
        {
          group = group + 1;
          if ( order * ( x1 - x0 ) <= 0.0 )
          {
            x1 = x0 + order * eps3;
          }
        }
      }
/*
  Check for isolated root.
*/
      else
      {
        xu = 1.0;

        if ( p == q )
        {
          rv6[p] = 1.0;
          for ( k = 0; k < n; k++ )
          {
            z[k+r*n] = 0.0;
          }
          for ( k = p; k <= q; k++ )
          {
            z[k+r*n] = rv6[k] * xu;
          }
          x0 = x1;
          continue;
        }

        norm = fabs ( d[p] );

        for ( i = p + 1; i <= q; i++ )
        {
          norm = r8_max ( norm, fabs ( d[i] ) + fabs ( e[i] ) );
        }
/*
  EPS2 is the criterion for grouping,
  EPS3 replaces zero pivots and equal roots are modified by EPS3,
  EPS4 is taken very small to avoid overflow.
*/
        eps2 = 0.001 * norm;
        eps3 = fabs ( norm ) * r8_epsilon ( );
        uk = q - p + 1;
        eps4 = uk * eps3;
        uk = eps4 / sqrt ( uk );
        s = p;
        group = 0;
      }
/*
  Elimination with interchanges and initialization of vector.
*/
      v = 0.0;

      for ( i = p; i <= q; i++ )
      {
        rv6[i] = uk;

        if ( i == p )
        {
          u = d[i] - x1 - xu * v;
          if ( i != q )
          {
            v = e[i+1];
          }
        }
        else if ( fabs ( u ) <= fabs ( e[i] ) )
        {
          xu = u / e[i];
          rv4[i] = xu;
          rv1[i-1] = e[i];
          rv2[i-1] = d[i] - x1;
          rv3[i-1] = 0.0;
          if ( i != q )
          {
            rv3[i-1] = e[i+1];
          }
          u = v - xu * rv2[i-1];
          v = - xu * rv3[i-1];
        }
        else
        {
          xu = e[i] / u;
          rv4[i] = xu;
          rv1[i-1] = u;
          rv2[i-1] = v;
          rv3[i-1] = 0.0;

          u = d[i] - x1 - xu * v;
          if ( i != q )
          {
            v = e[i+1];
          }
        }
      }

      if ( u == 0.0 )
      {
        u = eps3;
      }

      rv1[q] = u;
      rv2[q] = 0.0;
      rv3[q] = 0.0;
/*
  Back substitution.
*/
      while ( true )
      {
        for ( i = q; p <= i; i-- )
        {
          rv6[i] = ( rv6[i] - u * rv2[i] - v * rv3[i] ) / rv1[i];
          v = u;
          u = rv6[i];
        }
/*
  Orthogonalize with respect to previous members of group.
*/
        j = r;

        for ( jj = 1; jj <= group; jj++ )
        {
          while ( true )
          {
            j = j - 1;

            if ( ind[j] == tag )
            {
              break;
            }
          }

          xu = 0.0;
          for ( k = p; k <= q; k++ )
          {
            xu = xu + rv6[k] * z[k+j*n];
          }
          for ( k = p; k <= q; k++ )
          {
            rv6[k] = rv6[k] - xu * z[k+j*n];
          }
        }

        norm = 0.0;
        for ( k = p; k <= q; k++ )
        {
          norm = norm + fabs ( rv6[k] );
        }
/*
  Normalize so that sum of squares is 1.
*/
        if ( 1.0 <= norm )
        {
          u = 0.0;
          for( i = p; i <= q; i++ )
          {
            u = pythag ( u, rv6[i] );
          }
          xu = 1.0 / u;

          for ( k = 0; k < n; k++ )
          {
            z[k+r*n] = 0.0;
          }
          for ( k = p; k <= q; k++ )
          {
            z[k+r*n] = rv6[k] * xu;
          }
          x0 = x1;
          break;
        }
/*
  Set error: non-converged eigenvector.
*/
        else if ( 5 <= its )
        {
          ierr = - r;
          xu = 0.0;
          for ( k = 0; k < n; k++ )
          {
            z[k+r*n] = 0.0;
          }
          for ( k = p; k <= q; k++ )
          {
            z[k+r*n] = rv6[k] * xu;
          }
          x0 = x1;
          break;
        }
        else
        {
          if ( norm == 0.0 )
          {
            rv6[s] = eps4;
            s = s + 1;
            if ( q < s )
            {
              s = p;
            }
          }
          else
          {
            xu = eps4 / norm;
            for ( k = p; k <= q; k++ )
            {
              rv6[k] = rv6[k] * xu;
            }
          }
/*
  If RV1(I-1) == E(I), a row interchange was performed earlier in the
  triangularization process.
*/
          for ( i = p + 1; i <= q; i++ )
          {
            u = rv6[i];

            if ( rv1[i-1] == e[i] )
            {
              u = rv6[i-1];
              rv6[i-1] = rv6[i];
            }
            rv6[i] = u - rv4[i] * rv6[i-1];
          }
          its = its + 1;
        }
      }
    }

    if ( n <= q )
    {
      break;
    }
  }

  free ( rv1 );
  free ( rv2 );
  free ( rv3 );
  free ( rv4 );
  free ( rv6 );

  return ierr;
}
/******************************************************************************/

int tql1 ( int n, double d[], double e[] )

/******************************************************************************/
/*
  Purpose:

    TQL1 computes all eigenvalues of a real symmetric tridiagonal matrix.

  Discussion:

    TQL1 finds the eigenvalues of a symmetric tridiagonal
    matrix by the QL method.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  References:

    Bowdler, Martin, Reinsch, James Wilkinson,
    Numerische Mathematik,
    Volume 11, 1968, pages 293-306.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, is the order of the matrix.

    Input/output, double D[N].
    On input, the diagonal elements of the matrix.
    On output, the eigenvalues in ascending order.
    If an error exit is made, the eigenvalues are correct and
    ordered for indices 1, 2,... IERR-1, but may not be
    the smallest eigenvalues.

    Input/output, double E[N].  On input, E(2:N) contains the
    subdiagonal elements of the input matrix, and E(1) is arbitrary.
    On output, E has been destroyed.

    Output, int TQL1, error flag.
    0, normal return,
    J, if the J-th eigenvalue has not been determined after 30 iterations.
*/
{
  double c;
  double c2;
  double c3 = 0.0 ;
  double dl1;
  double el1;
  double f;
  double g;
  double h;
  int i;
  int ierr;
  int j;
  int l;
  int l1;
  int l2;
  int m;
  double p;
  double r;
  double s;
  double s2 = 0.0 ;
  double tst1;
  double tst2;

  ierr = 0;

  if ( n == 1 )
  {
    return ierr;
  }

  for ( i = 1; i < n; i++ )
  {
    e[i-1] = e[i];
  }

  f = 0.0;
  tst1 = 0.0;
  e[n-1] = 0.0;

  for ( l = 0; l < n; l++ )
  {
    j = 0;
    h = fabs ( d[l] ) + fabs ( e[l] );
    tst1 = r8_max ( tst1, h );
/*
  Look for a small sub-diagonal element.
*/
    for ( m = l; m < n; m++ )
    {
      tst2 = tst1 + fabs ( e[m] );

      if ( tst2 == tst1 )
      {
        break;
      }
    }

    if ( m != l )
    {
      while ( true )
      {
        if ( 30 <= j )
        {
          ierr = l + 1;
          return ierr;
        }

        j = j + 1;
/*
  Form the shift.
*/
        l1 = l + 1;
        l2 = l1 + 1;
        g = d[l];
        p = ( d[l1] - g ) / ( 2.0 * e[l] );
        r = pythag ( p, 1.0 );
        d[l] = e[l] / ( p + fabs ( r ) * r8_sign ( p ) );
        d[l1] = e[l] * ( p + fabs ( r ) * r8_sign ( p ) );
        dl1 = d[l1];
        h = g - d[l];

        for ( i = l2; i < n; i++ )
        {
          d[i] = d[i] - h;
        }

        f = f + h;
/*
  QL transformation.
*/
        p = d[m];
        c = 1.0;
        c2 = c;
        el1 = e[l1];
        s = 0.0;

        for ( i = m - 1; l <= i; i-- )
        {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = pythag ( p, e[i] );
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * ( c * g + s * d[i] );
        }

        p = - s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
        tst2 = tst1 + fabs ( e[l] );

        if ( tst2 <= tst1 )
        {
          break;
        }

      }
    }

    p = d[l] + f;
/*
  Order the eigenvalues.
*/
    for ( i = l; 0 <= i; i-- )
    {
      if ( i == 0 )
      {
        d[i] = p;
      }
      else if ( d[i-1] <= p )
      {
        d[i] = p;
        break;
      }
      else
      {
        d[i] = d[i-1];
      }
    }
  }

  return ierr;
}
/******************************************************************************/

int tql2 ( int n, double d[], double e[], double z[] )

/******************************************************************************/
/*
  Purpose:

    TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.

  Discussion:

    TQL2 finds the eigenvalues and eigenvectors of a symmetric
    tridiagonal matrix by the QL method.  The eigenvectors of a full
    symmetric matrix can also be found if TRED2 has been used to reduce this
    full matrix to tridiagonal form.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Bowdler, Martin, Reinsch, Wilkinson,
    TQL2,
    Numerische Mathematik,
    Volume 11, pages 293-306, 1968.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D[N].  On input, the diagonal elements of
    the matrix.  On output, the eigenvalues in ascending order.  If an error
    exit is made, the eigenvalues are correct but unordered for indices
    1,2,...,IERR-1.

    Input/output, double E[N].  On input, E(2:N) contains the
    subdiagonal elements of the input matrix, and E(1) is arbitrary.
    On output, E has been destroyed.

    Input, double Z[N*N].  On input, the transformation matrix
    produced in the reduction by TRED2, if performed.  If the eigenvectors of
    the tridiagonal matrix are desired, Z must contain the identity matrix.
    On output, Z contains the orthonormal eigenvectors of the symmetric
    tridiagonal (or full) matrix.  If an error exit is made, Z contains
    the eigenvectors associated with the stored eigenvalues.

    Output, int TQL2, error flag.
    0, normal return,
    J, if the J-th eigenvalue has not been determined after
    30 iterations.
*/
{
  double c;
  double c2;
  double c3 = 0.0;
  double dl1;
  double el1;
  double f;
  double g;
  double h;
  int i;
  int ierr;
  int ii;
  int j;
  int k;
  int l;
  int l1;
  int l2;
  int m;
  int mml;
  double p;
  double r;
  double s;
  double s2 = 0.0 ;
  double t;
  double tst1;
  double tst2;

  ierr = 0;

  if ( n == 1 )
  {
    return ierr;
  }

  for ( i = 1; i < n; i++ )
  {
    e[i-1] = e[i];
  }

  f = 0.0;
  tst1 = 0.0;
  e[n-1] = 0.0;

  for ( l = 0; l < n; l++ )
  {
    j = 0;
    h = fabs ( d[l] ) + fabs ( e[l] );
    tst1 = r8_max ( tst1, h );
/*
  Look for a small sub-diagonal element.
*/
    for ( m = l; m < n; m++ )
    {
      tst2 = tst1 + fabs ( e[m] );
      if ( tst2 == tst1 )
      {
        break;
      }
    }

    if ( m != l )
    {
      for ( ; ; )
      {
        if ( 30 <= j )
        {
          ierr = l + 1;
          return ierr;
        }

        j = j + 1;
/*
  Form shift.
*/
        l1 = l + 1;
        l2 = l1 + 1;
        g = d[l];
        p = ( d[l1] - g ) / ( 2.0 * e[l] );
        r = pythag ( p, 1.0 );
        d[l] = e[l] / ( p + r8_sign ( p ) * fabs ( r ) );
        d[l1] = e[l] * ( p + r8_sign ( p ) * fabs ( r ) );
        dl1 = d[l1];
        h = g - d[l];
        for ( i = l2; i < n; i++ )
        {
          d[i] = d[i] - h;
        }
        f = f + h;
/*
  QL transformation.
*/
        p = d[m];
        c = 1.0;
        c2 = c;
        el1 = e[l1];
        s = 0.0;
        mml = m - l;

        for ( ii = 1; ii <= mml; ii++ )
        {
          c3 = c2;
          c2 = c;
          s2 = s;
          i = m - ii;
          g = c * e[i];
          h = c * p;
          r = pythag ( p, e[i] );
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * ( c * g + s * d[i] );
/*
  Form vector.
*/
          for ( k = 0; k < n; k++ )
          {
            h = z[k+(i+1)*n];
            z[k+(i+1)*n] = s * z[k+i*n] + c * h;
            z[k+i*n] = c * z[k+i*n] - s * h;
          }
        }
        p = - s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
        tst2 = tst1 + fabs ( e[l] );

        if ( tst2 <= tst1 )
        {
          break;
        }
      }
    }
    d[l] = d[l] + f;
  }
/*
  Order eigenvalues and eigenvectors.
*/
  for ( ii = 1; ii < n; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i];
    for ( j = ii; j < n; j++ )
    {
      if ( d[j] < p )
      {
        k = j;
        p = d[j];
      }
    }

    if ( k != i )
    {
      d[k] = d[i];
      d[i] = p;
      for ( j = 0; j < n; j++ )
      {
        t        = z[j+i*n];
        z[j+i*n] = z[j+k*n];
        z[j+k*n] = t;
      }
    }
  }
  return ierr;
}
/******************************************************************************/

int tqlrat ( int n, double d[], double e2[] )

/******************************************************************************/
/*
  Purpose:

    TQLRAT computes all eigenvalues of a real symmetric tridiagonal matrix.

  Discussion:

    TQLRAT finds the eigenvalues of a symmetric
    tridiagonal matrix by the rational QL method.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Christian Reinsch,
    Algorithm 464, TQLRAT,
    Communications of the ACM,
    Volume 16, page 689, 1973.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D[N].  On input, D contains the diagonal
    elements of the matrix.  On output, D contains the eigenvalues in ascending
    order.  If an error exit was made, then the eigenvalues are correct
    in positions 1 through IERR-1, but may not be the smallest eigenvalues.

    Input/output, double E2[N], contains in positions 2 through N 
    the squares of the subdiagonal elements of the matrix.  E2(1) is
    arbitrary.  On output, E2 has been overwritten by workspace
    information.

    Output, int TQLRAT, error flag.
    0, for no error,
    J, if the J-th eigenvalue could not be determined after 30 iterations.
*/
{
  double b = 0.0 ;
  double c = 0.0 ;
  double f;
  double g;
  double h;
  int i;
  int ierr;
  int ii;
  int j;
  int l;
  int l1;
  int m;
  int mml;
  double p;
  double r;
  double s;
  double t;

  ierr = 0;

  if ( n == 1 )
  {
    return ierr;
  }

  for ( i = 1; i < n; i++ )
  {
    e2[i-1] = e2[i];
  }

  f = 0.0;
  t = 0.0;
  e2[n-1] = 0.0;

  for ( l = 0; l < n; l++ )
  {
     j = 0;
     h = fabs ( d[l] ) + sqrt ( e2[l] );

     if ( t <= h )
     {
       t = h;
       b = fabs ( t ) * r8_epsilon ( );
       c = b * b;
     }
/*
  Look for small squared sub-diagonal element.
*/
    for ( m = l; m < n; m++ )
    {  
      if ( e2[m] <= c )
      {
        break;
      }
    }

    if ( m != l )
    {
      for ( ; ; )
      {
        if ( 30 <= j )
        {
          ierr = l + 1;
          return ierr;
        }

        j = j + 1;
/*
  Form shift.
*/
        l1 = l + 1;
        s = sqrt ( e2[l] );
        g = d[l];
        p = ( d[l1] - g ) / ( 2.0 * s );
        r = pythag ( p, 1.0 );
        d[l] = s / ( p + fabs ( r ) * r8_sign ( p ) );
        h = g - d[l];
        for ( i = l1; i < n; i++ )
        {
          d[i] = d[i] - h;
        }
        f = f + h;
/*
  Rational QL transformation.
*/
        g = d[m];
        if ( g == 0.0 )
        {
          g = b;
        }

        h = g;
        s = 0.0;
        mml = m - l;

        for ( ii = 1; ii <= mml; ii++ )
        {
          i = m - ii;
          p = g * h;
          r = p + e2[i];
          e2[i+1] = s * r;
          s = e2[i] / r;
          d[i+1] = h + s * ( h + d[i] );
          g = d[i] - e2[i] / g;
          if ( g == 0.0 )
          {
            g = b;
          }
          h = g * p / r;
        }
        e2[l] = s * g;
        d[l] = h;
/*
  Guard against underflow in convergence test.
*/
        if ( h == 0.0 )
        {
          break;
        }

        if ( fabs ( e2[l] ) <= fabs ( c / h ) )
        {
          break;
        }

        e2[l] = h * e2[l];

        if ( e2[l] == 0.0 )
        {
          break;
        }
      }
    }

    p = d[l] + f;
/*
  Order the eigenvalues.
*/
    for ( i = l; 0 <= i; i-- )
    {
      if ( i == 0 )
      {
        d[i] = p;
        break;
      }
      else if ( d[i-1] <= p )
      {
        d[i] = p;
        break;
      }
      d[i] = d[i-1];
    }
  }

  return ierr;
}
/******************************************************************************/

void trbak1 ( int n, double a[], double e[], int m, double z[] )

/******************************************************************************/
/*
  Purpose:

    TRBAK1 determines eigenvectors by undoing the TRED1 transformation.

  Discussion:

    TRBAK1 forms the eigenvectors of a real symmetric
    matrix by back transforming those of the corresponding
    symmetric tridiagonal matrix determined by TRED1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N,N], contains information about the orthogonal
    transformations used in the reduction by TRED1 in its strict lower
    triangle.

    Input, double E[N], the subdiagonal elements of the tridiagonal
    matrix in E(2:N).  E(1) is arbitrary.

    Input, int M, the number of eigenvectors to be back
    transformed.

    Input/output, double Z[N,M].  On input, the eigenvectors to be
    back transformed.  On output, the transformed eigenvectors.
*/
{
  int i;
  int j;
  int k;
  int l;
  double s;

  if ( m <= 0 )
  {
    return;
  }

  if ( n <= 1 )
  {
    return;
  }

  for ( i = 1; i < n; i++ )
  {
    l = i - 1;

    if ( e[i] != 0.0 )
    {
      for ( j = 0; j < m; j++ )
      {
        s = 0.0;
        for ( k = 0; k < i; k++ )
        {
          s = s + a[i+k*n] * z[k+j*n];
        }
        s = ( s / a[i+l*n] ) / e[i];
        for ( k = 0; k < i; k++ )
        {
          z[k+j*n] = z[k+j*n] + s * a[i+k*n];
        }
      }
    }
  }

  return;
}
/******************************************************************************/

void trbak3 ( int n, int nv, double a[], int m, double z[] )

/******************************************************************************/
/*
  Purpose:

    TRBAK3 determines eigenvectors by undoing the TRED3 transformation.

  Discussion:

    TRBAK3 forms the eigenvectors of a real symmetric
    matrix by back transforming those of the corresponding
    symmetric tridiagonal matrix determined by TRED3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

   Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int NV, the dimension of the array paramater A,
    which must be at least N*(N+1)/2.

    Input, double A[NV], information about the orthogonal
    transformations used in the reduction by TRED3.

    Input, int M, the number of eigenvectors to be back
    transformed.

    Input/output, double Z[N,M].  On input, the eigenvectors to be 
    back transformed.  On output, the transformed eigenvectors.
*/
{
  double h;
  int i;
  int ik;
  int iz;
  int j;
  int k;
  double s;

  if ( m == 0 )
  {
    return;
  }

  for ( i = 1; i < n; i++ )
  {
    iz = ( i * ( i + 1 ) ) / 2;
    ik = iz + i;
    h = a[ik];

    if ( h != 0.0 )
    {
      for ( j = 0; j < m; j++ )
      {
        s = 0.0;
        ik = iz - 1;

        for ( k = 0; k < i; k++ )
        {
          ik = ik + 1;
          s = s + a[ik] * z[k+j*n];
        }

        s = ( s / h ) / h;
        ik = iz - 1;

        for ( k = 0; k < i; k++ )
        {
          ik = ik + 1;
          z[k+j*n] = z[k+j*n] - s * a[ik];
        }
      }
    }
  }

  return;
}
/******************************************************************************/

void tred1 ( int n, double a[], double d[], double e[], double e2[] )

/******************************************************************************/
/*
  Purpose:

    TRED1 transforms a real symmetric matrix to symmetric tridiagonal form.

  Discussion:

    TRED1 reduces a real symmetric matrix to a symmetric
    tridiagonal matrix using orthogonal similarity transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Martin, Reinsch, Wilkinson,
    TRED1,
    Numerische Mathematik,
    Volume 11, pages 181-195, 1968.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix A.

    Input/output, double A[N*N], on input, contains the real
    symmetric matrix.  Only the lower triangle of the matrix need be supplied.
    On output, A contains information about the orthogonal transformations
    used in the reduction in its strict lower triangle.
    The full upper triangle of A is unaltered.

    Output, double D[N], contains the diagonal elements of the
    tridiagonal matrix.

    Output, double E[N], contains the subdiagonal elements of the
    tridiagonal matrix in its last N-1 positions.  E(1) is set to zero.

    Output, double E2[N], contains the squares of the corresponding
    elements of E.  E2 may coincide with E if the squares are not needed.
*/
{
  double f;
  double g;
  double h;
  int i;
  int j;
  int k;
  int l;
  double scale;

  for ( j = 0; j < n; j++ )
  {
    d[j] = a[n-1+j*n];
  }

  for ( i = 0; i < n; i++ )
  {
    a[n-1+i*n] = a[i+i*n];
  }

  for ( i = n - 1; 0 <= i; i-- )
  {
    l = i - 1;
    h = 0.0;
/*
  Scale row.
*/
    scale = 0.0;
    for ( k = 0; k <= l; k++ )
    {
      scale = scale + fabs ( d[k] );
    }

    if ( scale == 0.0 )
    {
      for ( j = 0; j <= l; j++ )
      {
        d[j]     = a[l+j*n];
        a[l+j*n] = a[i+j*n];
        a[i+j*n] = 0.0;
      }

      e[i] = 0.0;
      e2[i] = 0.0;
      continue;
    }

    for ( k = 0; k <= l; k++ )
    {
      d[k] = d[k] / scale;
    }

    for ( k = 0; k <= l; k++ )
    {
      h = h + d[k] * d[k];
    }

    e2[i] = h * scale * scale;
    f = d[l];
    g = - sqrt ( h ) * r8_sign ( f );
    e[i] = scale * g;
    h = h - f * g;
    d[l] = f - g;

    if ( 0 <= l )
    {
/*
  Form A * U.
*/
      for ( k = 0; k <= l; k++ )
      {
        e[k] = 0.0;
      }

      for ( j = 0; j <= l; j++ )
      {
        f = d[j];
        g = e[j] + a[j+j*n] * f;

        for ( k = j + 1; k <= l; k++ )
        {
          g = g + a[k+j*n] * d[k];
          e[k] = e[k] + a[k+j*n] * f;
        }
        e[j] = g;
      }
/*
  Form P.
*/
      f = 0.0;
      for ( j = 0; j <= l; j++ )
      {
        e[j] = e[j] / h;
        f = f + e[j] * d[j];
      }

      h = f / ( h + h );
/*
  Form Q.
*/
      for ( j = 0; j <= l; j++ )
      {
        e[j] = e[j] - h * d[j];
      }
/*
  Form reduced A.
*/
      for ( j = 0; j <= l; j++ )
      {
        f = d[j];
        g = e[j];
        for ( k = j; k <= l; k++ )
        {
          a[k+j*n] = a[k+j*n] - f * e[k] - g * d[k];
        }
      }
    }

    for ( j = 0; j <= l; j++ )
    {
      f        = d[j];
      d[j]     = a[l+j*n];
      a[l+j*n] = a[i+j*n];
      a[i+j*n] = f * scale;
    }
  }
  return;
}
/******************************************************************************/

void tred2 ( int n, double a[], double d[], double e[], double z[] )

/******************************************************************************/
/*
  Purpose:

    TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.

  Discussion:

    TRED2 reduces a real symmetric matrix to a
    symmetric tridiagonal matrix using and accumulating
    orthogonal similarity transformations.

    A and Z may coincide, in which case a single storage area is used
    for the input of A and the output of Z.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 November 2012

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Martin, Reinsch, Wilkinson,
    TRED2,
    Numerische Mathematik,
    Volume 11, pages 181-195, 1968.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the real symmetric input matrix.  Only the
    lower triangle of the matrix need be supplied.

    Output, double D[N], the diagonal elements of the tridiagonal
    matrix.

    Output, double E[N], contains the subdiagonal elements of the
    tridiagonal matrix in E(2:N).  E(1) is set to zero.

    Output, double Z[N*N], the orthogonal transformation matrix
    produced in the reduction.
*/
{
  double f;
  double g;
  double h;
  double hh;
  int i;
  int j;
  int k;
  int l;
  double scale;

  for ( j = 0; j < n; j++ )
  {
    for ( i = j; i < n; i++ )
    {
      z[i+j*n] = a[i+j*n];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    d[j] = a[n-1+j*n];
  }

  for ( i = n - 1; 1 <= i; i-- )
  {
    l = i - 1;
    h = 0.0;
/*
  Scale row.
*/
    scale = 0.0;
    for ( k = 0; k <= l; k++ )
    {
      scale = scale + fabs ( d[k] );
    }

    if ( scale == 0.0 )
    {
      e[i] = d[l];

      for ( j = 0; j <= l; j++ )
      {
        d[j]     = z[l+j*n];
        z[i+j*n] = 0.0;
        z[j+i*n] = 0.0;
      }
      d[i] = 0.0;
      continue;
    }

    for ( k = 0; k <= l; k++ )
    {
      d[k] = d[k] / scale;
    }

    h = 0.0;
    for ( k = 0; k <= l; k++ )
    {
      h = h + d[k] * d[k];
    }

    f = d[l];
    g = - sqrt ( h ) * r8_sign ( f );
    e[i] = scale * g;
    h = h - f * g;
    d[l] = f - g;
/*
  Form A*U.
*/
    for ( k = 0; k <= l; k++ )
    {
      e[k] = 0.0;
    }

    for ( j = 0; j <= l; j++ )
    {
      f = d[j];
      z[j+i*n] = f;
      g = e[j] + z[j+j*n] * f;

      for ( k = j + 1; k <= l; k++ )
      {
        g = g + z[k+j*n] * d[k];
        e[k] = e[k] + z[k+j*n] * f;
      }
      e[j] = g;
    }
/*
  Form P.
*/
    for ( k = 0; k <= l; k++ )
    {
      e[k] = e[k] / h;
    }
    f = 0.0;
    for ( k = 0; k <= l; k++ )
    {
      f = f + e[k] * d[k];
    }
    hh = 0.5 * f / h;
/*
  Form Q.
*/
    for ( k = 0; k <= l; k++ )
    {
      e[k] = e[k] - hh * d[k];
    }
/*
  Form reduced A.
*/
    for ( j = 0; j <= l; j++ )
    {
      f = d[j];
      g = e[j];

      for ( k = j; k <= l; k++ )
      {
        z[k+j*n] = z[k+j*n] - f * e[k] - g * d[k];
      }
      d[j] = z[l+j*n];
      z[i+j*n] = 0.0;
    }
    d[i] = h;
  }
/*
  Accumulation of transformation matrices.
*/
  for ( i = 1; i < n; i++ )
  {
    l = i - 1;
    z[n-1+l*n] = z[l+l*n];
    z[l+l*n] = 1.0;
    h = d[i];

    if ( h != 0.0 )
    {
      for ( k = 0; k <= l; k++ )
      {
        d[k] = z[k+i*n] / h;
      }
      for ( j = 0; j <= l; j++ )
      {
        g = 0.0;
        for ( k = 0; k <= l; k++ )
        {
          g = g + z[k+i*n] * z[k+j*n];
        }
        for ( k = 0; k <= l; k++ )
        {
          z[k+j*n] = z[k+j*n] - g * d[k];
        }
      }
    }
    for ( k = 0; k <= l; k++ )
    {
      z[k+i*n] = 0.0;
    }
  }

  for ( j = 0; j < n; j++ )
  {
    d[j] = z[n-1+j*n];
  }

  for ( j = 0; j < n - 1; j++ )
  {
    z[n-1+j*n] = 0.0;
  }
  z[n-1+(n-1)*n] = 1.0;

  e[0] = 0.0;

  return;
}
/******************************************************************************/

void tred3 ( int n, int nv, double a[], double d[], double e[], double e2[] )

/******************************************************************************/
/*
  Purpose:

    TRED3: transform real symmetric packed matrix to symmetric tridiagonal form.

  Discussion:

    TRED3 reduces a real symmetric matrix, stored as
    a one-dimensional array, to a symmetric tridiagonal matrix
    using orthogonal similarity transformations.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 January 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt.

  Reference:

    Martin, Reinsch, James Wilkinson,
    TRED3,
    Numerische Mathematik,
    Volume 11, pages 181-195, 1968.

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input, int NV, the dimension of A, which must be at least
    (N*(N+1))/2.

    Input/output, double A[NV].  On input, the lower triangle of
    the real symmetric matrix, stored row-wise.  On output, information about
    the orthogonal transformations used in the reduction.

    Output, double D[N], the diagonal elements of the tridiagonal
    matrix.

    Output, double E[N], the subdiagonal elements of the tridiagonal
    matrix in E(2:N).  E(1) is set to zero.

    Output, double E2[N],  the squares of the corresponding
    elements of E.  E2 may coincide with E if the squares are not needed.
*/
{
  double f;
  double g;
  double h;
  double hh;
  int i;
  int iz;
  int j;
  int jk;
  int k;
  double scale;

  for ( i = n - 1; 0 <= i; i-- )
  {
    iz = ( i * ( i + 1 ) ) / 2 - 1;
    h = 0.0;
    scale = 0.0;
/*
  Scale row.
*/
    for ( k = 0; k < i; k++ )
    {
      iz = iz + 1;
      d[k] = a[iz];
      scale = scale + fabs ( d[k] );
    }

    if ( scale == 0.0 )
    {
      e[i] = 0.0;
      e2[i] = 0.0;
      d[i] = a[iz+1];
      a[iz+1] = scale * sqrt ( h );
      continue;
    }

    for ( k = 0; k < i; k++ )
    {
      d[k] = d[k] / scale;
      h = h + d[k] * d[k];
    }

    e2[i] = scale * scale * h;
    f = d[i-1];
    g = - sqrt ( h ) * r8_sign ( f );
    e[i] = scale * g;
    h = h - f * g;
    d[i-1] = f - g;
    a[iz] = scale * d[i-1];
    if ( i == 1 )
    {
      d[i] = a[iz+1];
      a[iz+1] = scale * sqrt ( h );
      continue;
    }

    jk = 0;

    for ( j = 0; j < i; j++ )
    {
      f = d[j];
      g = 0.0;

      for ( k = 0; k < j; k++ )
      {
        g = g + a[jk] * d[k];
        e[k] = e[k] + a[jk] * f;
        jk = jk + 1;
      }

      e[j] = g + a[jk] * f;
      jk = jk + 1;
    }
/*
  Form P.
*/
    for ( j = 0; j < i; j++ )
    {
      e[j] = e[j] / h;
    }
    f = 0.0;
    for ( j = 0; j < i; j++ )
    {
      f = f + e[j] * d[j];
    }
    hh = f / ( h + h );
/*
  Form Q.
*/
    for ( j = 0; j < i; j++ )
    {
      e[j] = e[j] - hh * d[j];
    }
    jk = 0;
/*
  Form reduced A.
*/
    for ( j = 0; j < i; j++ )
    {
      f = d[j];
      g = e[j];
      for ( k = 0; k <= j; k++ )
      {
        a[jk] = a[jk] - f * e[k] - g * d[k];
        jk = jk + 1;
      }
    }

    d[i] = a[iz+1];
    a[iz+1] = scale * sqrt ( h );
  }

  return;
}
/******************************************************************************/

int tridib ( int n, double *eps1, double d[], double e[], double e2[], 
  double *lb, double *ub, int m11, int m, double w[], int ind[] )

/******************************************************************************/
/*
  Purpose:

    TRIDIB computes some eigenvalues of a real symmetric tridiagonal matrix.

  Discussion:

    TRIDIB finds those eigenvalues of a tridiagonal symmetric matrix between 
    specified boundary indices, using bisection.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double *EPS1.  On input, an absolute error
    tolerance for the computed eigenvalues.  It should be chosen commensurate
    with relative perturbations in the matrix elements of the order of the
    relative machine precision.  If the input EPS1 is non-positive, it
    is reset for each submatrix to a default value, namely, minus the
    product of the relative machine precision and the 1-norm of the submatrix.

    Input, double D[N], the diagonal elements of the input matrix.

    Input, double E[N], the subdiagonal elements of the input matrix
    in E[1:N-1].  E[0] is arbitrary.

    Input/output, double E2[N].  On input, the squares of the
    corresponding elements of E.  E2(1) is arbitrary.  On output, elements of
    E2 corresponding to elements of E regarded as negligible, have been
    replaced by zero, causing the matrix to split into a direct sum of
    submatrices.  E2(1) is also set to zero.

    Output, double *LB, *UB, define an interval containing exactly
    the desired eigenvalues.

    Input, int M11, the lower boundary index for the desired eigenvalues.
    0 <= M11 < N.

    Input, int M, the number of eigenvalues desired.  The
    upper boundary index M22 is then obtained as M22 = M11 + M - 1.

    Output, double W[M], the eigenvalues between indices M11 and M22
    in ascending order.

    Output, int IND[M], the submatrix indices associated with
    the corresponding eigenvalues in W: 1 for eigenvalues belonging to the
    first submatrix from the top, 2 for those belonging to the second
    submatrix, and so on.

    Output, int TRIDIB, error flag.
    0, for normal return,
    3*N+1, if multiple eigenvalues at index M11 make unique selection
      impossible,
    3*N+2, if multiple eigenvalues at index M22 make unique selection
      impossible.
*/
{
  int i;
  int ierr;
  int isturm;
  int j;
  int k;
  int l;
  int m1;
  int m2;
  int m22;
  int p;
  int q;
  int r;
  double *rv4;
  double *rv5;
  int s;
  double t1;
  double t2;
  int tag;
  double tst1;
  double tst2;
  double u;
  double v;
  double x0;
  double x1;
  double xu;

  rv4 = ( double * ) malloc ( n * sizeof ( double ) );
  rv5 = ( double * ) malloc ( n * sizeof ( double ) );

  ierr = 0;
  tag = 0;
  xu = d[0];
  x0 = d[0];
  s = 0;
  u = 0.0;
/*
  Look for small sub-diagonal entries and determine an
  interval containing all the eigenvalues.
*/
  for ( i = 1; i <= n; i++ )
  {
    x1 = u;

    if ( i == n )
    {
      u = 0.0;
    }
    else
    {
      u = fabs ( e[i] );
    }

    xu = r8_min ( xu, d[i-1] - ( x1 + u ) );
    x0 = r8_max ( x0, d[i-1] + ( x1 + u ) );

    if ( 1 < i )
    {
      tst1 = fabs ( d[i-1] ) + fabs ( d[i-2] );
      tst2 = tst1 + fabs ( e[i-1] );
      if ( tst2 <= tst1 )
      {
        e2[i-1] = 0.0;
      }
    }
    else
    {
      e2[i-1] = 0.0;
    }
  }

  x1 = ( double ) n;
  x1 = x1 * r8_max ( fabs ( xu ), fabs ( x0 ) ) * r8_epsilon ( );
  xu = xu - x1;
  t1 = xu;
  x0 = x0 + x1;
  t2 = x0;
/*
  Determine an interval containing exactly the desired eigenvalues.
*/
  p = 1;
  q = n;
  m1 = m11 - 1;
  m22 = m1 + m;

  if ( m1 != 0 || m != n )
  {
    if ( m1 == 0 && m != n )
    {
       x0 = t2;
       isturm = 2;
    }
    else
    {
      isturm = 1;
    }

    while ( true )
    {
      v = x1;
      x1 = xu + ( x0 - xu ) * 0.5;

      if ( x1 == v )
      {
        ierr = 3 * n + isturm;
        free ( rv4 );
        free ( rv5 );
        *lb = t1;
        *ub = t2;
        return ierr;
      }

      s = sturm_sequence ( d, e, e2, n, p - 1, q - 1, x1 );

      if ( isturm == 1 )
      {
        if ( s < m1 )
        {
          xu = x1;
          continue;
        }
        else if ( m1 < s )
        {
          x0 = x1;
          continue;
        }
        xu = x1;
        t1 = x1;
        m22 = m1 + m;

        if ( m22 != n )
        {
          x0 = t2;
          isturm = 2;
          continue;
        }
      }
      else
      {
        if ( s < m22 )
        {
          xu = x1;
          continue;
        }
        else if ( m22 < s )
        {
          x0 = x1;
          continue;
        }
        t2 = x1;
      }
      break;
    }
  }

  q = 0;
  r = 0;
/*
  Establish and process next submatrix, refining interval by the
  Gerschgorin bounds.
*/
  while ( true )
  {
    if ( r == m )
    {
      free ( rv4 );
      free ( rv5 );
      *lb = t1;
      *ub = t2;
      return ierr;
    }

    tag = tag + 1;
    p = q + 1;
    xu = d[p-1];
    x0 = d[p-1];
    u = 0.0;

    for ( q = p; q <= n; q++ )
    {
      x1 = u;
      u = 0.0;
      v = 0.0;

      if ( q < n )
      {
        u = fabs ( e[q] );
        v = e2[q];
      }

      xu = r8_min ( d[q-1] - ( x1 + u ), xu );
      x0 = r8_max ( d[q-1] + ( x1 + u ), x0 );

      if ( v == 0.0 )
      {
        break;
      }
    }

    x1 = r8_max ( fabs ( xu ), fabs ( x0 ) ) * r8_epsilon ( );

    if ( *eps1 <= 0.0 )
    {
      *eps1 = - x1;
    }
/*
  Check for isolated root within interval.
*/
    if ( p == q )
    {
      if ( d[p-1] < t1 || t2 <= d[p-1] )
      {
        if ( n <= q ) 
        {
          break;
        }
        continue;
      }
      m1 = p;
      m2 = p;
      rv5[p-1] = d[p-1];
    }
    else
    {
      x1 = x1 * ( q - p + 1 );
      *lb = r8_max ( t1, xu - x1 );
      *ub = r8_min ( t2, x0 + x1 );

      x1 = *lb;
      s = sturm_sequence ( d, e, e2, n, p - 1, q - 1, x1 );
      m1 = s + 1;

      x1 = *ub;
      s = sturm_sequence ( d, e, e2, n, p - 1, q - 1, x1 );
      m2 = s;

      if ( m2 < m1 )
      {
        if ( n <= q ) 
        {
          break;
        }
        continue;
      }
/*
  Find roots by bisection.
*/
      x0 = *ub;
      isturm = 5;

      for ( i = m1; i <= m2; i++ )
      {
        rv5[i-1] = *ub;
        rv4[i-1] = *lb;
      }
/*
  Loop for the K-th eigenvalue.
*/
      k = m2;

      while ( true )
      {
        xu = *lb;

        for ( i = k; m1 <= i; i-- )
        {
          if ( xu < rv4[i-1] )
          {
            xu = rv4[i-1];
            break;
          }
        }
        x0 = r8_min ( x0, rv5[k-1] );
/*
  Next bisection step.
*/
        while ( true )
        {
          x1 = ( xu + x0 ) * 0.5;
  
          if ( ( x0 - xu ) <= fabs ( *eps1 ) )
          {
            break;
          }
          tst1 = 2.0 * ( fabs ( xu ) + fabs ( x0 ) );
          tst2 = tst1 + ( x0 - xu );

          if ( tst2 == tst1 )
          {
            break;
          }
          s = sturm_sequence ( d, e, e2, n, p - 1, q - 1, x1 );
/*
  Refine intervals.
*/
          if ( k <= s )
          {
            x0 = x1;
          }
          else
          {
            xu = x1;
            if ( m1 <= s )
            {
              rv4[s] = x1;
              rv5[s-1] = r8_min ( rv5[s-1], x1 );
            }
            else
            {
              rv4[m1-1] = x1;
            }
          }
        }
/*
  K-th eigenvalue found.
*/
        rv5[k-1] = x1;
        k = k - 1;

        if ( k < m1 )
        {
          break;
        }
      }
    }
/*
  Order eigenvalues tagged with their submatrix associations.
*/
    s = r;
    r = r + m2 - m1 + 1;
    j = 1;
    k = m1;

    for ( l = 1; l <= r; l++ )
    {
      if ( j <= s )
      { 
        if ( m2 < k )
        {
          break;
        }

        if ( w[l-1] <= rv5[k-1] )
        {
          j = j + 1;
          continue;
        }
        for ( i = l + s - j; l <= i; i-- )
        {
          w[i] = w[i-1];
          ind[i] = ind[i-1];
        }
      }
      w[l-1] = rv5[k-1];
      ind[l-1] = tag;
      k = k + 1;
    }

    if ( n <= q ) 
    {
      break;
    }
  }
  free ( rv4 );
  free ( rv5 );
  *lb = t1;
  *ub = t2;
  return ierr;
}
/******************************************************************************/

int tsturm ( int n, double *eps1, double d[], double e[], double e2[], 
  double lbb, double ubb, int mm, int *m, double w[], double z[] )

/******************************************************************************/
/*
  Purpose:

    TSTURM computes some eigenvalues/vectors, real symmetric tridiagonal matrix.

  Discussion:

    TSTURM finds those eigenvalues of a tridiagonal symmetric matrix which 
    lie in a specified interval and their associated eigenvectors, using 
    bisection and inverse iteration.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 March 2018

  Author:

    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
    Klema, Moler.
    C version by John Burkardt

  Reference:

    James Wilkinson, Christian Reinsch,
    Handbook for Automatic Computation,
    Volume II, Linear Algebra, Part 2,
    Springer, 1971,
    ISBN: 0387054146,
    LC: QA251.W67.

    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
    Matrix Eigensystem Routines, EISPACK Guide,
    Lecture Notes in Computer Science, Volume 6,
    Springer Verlag, 1976,
    ISBN13: 978-3540075462,
    LC: QA193.M37.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double *EPS1.  On input, an absolute error
    tolerance for the computed eigenvalues.  It should be chosen commensurate
    with relative perturbations in the matrix elements of the order of the
    relative machine precision.  If the input EPS1 is non-positive, it
    is reset for each submatrix to a default value, namely, minus the
    product of the relative machine precision and the 1-norm of the submatrix.

    Input, double D[N], the diagonal elements of the input matrix.

    Input, double E[N], the subdiagonal elements of the input matrix
    in E(2:N).  E(1) is arbitrary.

    Input/output, double E2(N).  On input, the squares of the
    corresponding elements of E.  E2(1) is arbitrary.  On output, elements of
    E2 corresponding to elements of E regarded as negligible have been
    replaced by zero, causing the matrix to split into a direct sum of
    submatrices.  E2(1) is also set to zero.

    Input, double LBB, UBB, define the interval to be searched for
    eigenvalues.  If LB is not less than UB, no eigenvalues will be found.

    Input, int MM, an upper bound for the number of
    eigenvalues in the interval.  If more than MM eigenvalues are determined
    to lie in the interval, an error return is made with no values or vectors
    found.

    Output, int *M, the number of eigenvalues determined to lie
    in (LB, UB).

    Output, double W[M], the eigenvalues in ascending order if the
    matrix does not split.  If the matrix splits, the eigenvalues are in
    ascending order for each submatrix.  If a vector error break; is made, W
    contains those values already found.

    Output, double Z[N,MM], the associated set of orthonormal
    eigenvectors.  If an error break; is made, Z contains those vectors already
    found.

    Output, int TSTURM, error flag.
    0, normal return.
    3*N+1, if M exceeds MM.
    4*N+R, if the eigenvector corresponding to the R-th
      eigenvalue fails to converge in 5 iterations.
*/
{
  double eps2;
  double eps3;
  double eps4;
  int group;
  int i;
  int ierr;
  int its;
  int j;
  int k;
  double lb;
  int m1;
  int m2;
  double norm;
  int p;
  int q;
  int r;
  double *rv1;
  double *rv2;
  double *rv3;
  double *rv4;
  double *rv5;
  double *rv6;
  int s;
  double tst1;
  double tst2;
  double u;
  double ub;
  double uk;
  double v;
  double x0;
  double x1;
  double xu;

  ierr = 0;
  s = 0;
  lb = lbb;
  ub = ubb;
/*
  Look for small sub-diagonal entries.
*/
  e2[0] = 0.0;
  for ( i = 1; i < n; i++ )
  {
    tst1 = fabs ( d[i] ) + fabs ( d[i-1] );
    tst2 = tst1 + fabs ( e[i] );
    if ( tst2 <= tst1 )
    {
      e2[i] = 0.0;
    }
  }
/*
  Determine the number of eigenvalues in the interval.
*/
  p = 0;
  q = n - 1;
  x1 = ub;
  s = sturm_sequence ( d, e, e2, n, p, q, x1 );
  *m = s;

  x1 = lb;
  s = sturm_sequence ( d, e, e2, n, p, q, x1 );
  *m = *m - s;

  if ( mm < *m )
  {
    ierr = 3 * n + 1;
    return ierr;
  }
/*
  Allocate work arrays.
*/
  rv1 = ( double * ) malloc ( n * sizeof ( double ) );
  rv2 = ( double * ) malloc ( n * sizeof ( double ) );
  rv3 = ( double * ) malloc ( n * sizeof ( double ) );
  rv4 = ( double * ) malloc ( n * sizeof ( double ) );
  rv5 = ( double * ) malloc ( n * sizeof ( double ) );
  rv6 = ( double * ) malloc ( n * sizeof ( double ) );

  q = -1;
  r = 0;
/*
  Establish and process next submatrix, refining interval by the
  Gerschgorin bounds.
*/
  while ( q < n - 1 )
  {
    if ( r == *m )
    {
      free ( rv1 );
      free ( rv2 );
      free ( rv3 );
      free ( rv4 );
      free ( rv5 );
      free ( rv6 );
      return ierr;
    }

    p = q + 1;
    xu = d[p];
    x0 = d[p];
    u = 0.0;

    for ( q = p; q < n; q++ )
    {
      x1 = u;
      u = 0.0;
      v = 0.0;

      if ( q < n - 1 )
      {
        u = fabs ( e[q+1] );
        v = e2[q+1];
      }

      xu = r8_min ( d[q] - ( x1 + u ), xu );
      x0 = r8_max ( d[q] + ( x1 + u ), x0 );

      if ( v == 0.0 )
      {
        break;
      }
    }

    x1 = r8_max ( fabs ( xu ), fabs ( x0 ) ) * r8_epsilon ( );

    if ( *eps1 <= 0.0 )
    {
      *eps1 = - x1;
    }
/*
  Check for isolated root within interval.
*/
    if ( p == q )
    {
      if ( d[p] < lbb || ubb <= d[p] )
      {
        continue;
      }

      r = r + 1;
      for ( i = 0; i < n; i++ )
      {
        z[i+(r-1)*n] = 0.0; 
      }
      w[r] = d[p];
      z[p+(r-1)*n] = 1.0;
      continue;
    }

    u = q - p + 1;
    x1 = u * x1;
    lb = r8_max ( lbb, xu - x1 );
    ub = r8_min ( ubb, x0 + x1 );
    x1 = lb;
    s = sturm_sequence ( d, e, e2, n, p, q, x1 );
    m1 = s;

    x1 = ub;
    s = sturm_sequence ( d, e, e2, n, p, q, x1 );
    m2 = s - 1;

    if ( m2 < m1 )
    {
      continue;
    }
/*
  Find roots by bisection.
*/
    x0 = ub;
    for ( i = m1; i <= m2; i++ )
    {
      rv5[i] = ub;
      rv4[i] = lb;
    }
/*
  Loop for K-th eigenvalue.
*/
    k = m2;

    while ( m1 <= k )
    {
      xu = lb;

      for ( i = k; m1 <= i; i-- )
      {
        if ( xu < rv4[i] )
        {
          xu = rv4[i];
          break;
        }
      }

      x0 = r8_min ( x0, rv5[k] );
/*
  Next bisection step.
*/
      while ( true )
      {
        x1 = ( xu + x0 ) * 0.5;

        if ( ( x0 - xu ) <= fabs ( *eps1 ) )
        {
          break;
        }

        tst1 = 2.0 * ( fabs ( xu ) + fabs ( x0 ) );
        tst2 = tst1 + ( x0 - xu );

        if ( tst2 == tst1 )
        {
          break;
        }

        s = sturm_sequence ( d, e, e2, n, p, q, x1 );
/*
  Refine intervals.
*/
        if ( k + 1 <= s )
        {
          x0 = x1;
          continue;
        }

        xu = x1;

        if ( m1 + 1 <= s )
        {
          rv4[s] = x1;
          if ( x1 < rv5[s-1] )
          {
            rv5[s-1] = x1;
          }
        }
        else
        {
          rv4[m1] = x1;
        }
      }
/*
  K-th eigenvalue found.
*/
      rv5[k] = x1;
      k = k - 1;
    }
/*
  Find vectors by inverse iteration.
*/
    norm = fabs ( d[p] );

    for ( i = p + 1; i <= q; i++ )
    {
      norm = r8_max ( norm, fabs ( d[i] ) + fabs ( e[i] ) );
    }
/*
  EPS2 is the criterion for grouping,
  EPS3 replaces zero pivots and equal roots are modified by eps3,
  EPS4 is taken very small to avoid overflow.
*/
    eps2 = 0.001 * norm;
    eps3 = fabs ( norm ) * r8_epsilon ( );
    uk = q - p + 1;
    eps4 = uk * eps3;
    uk = eps4 / sqrt ( uk );
    group = 0;
    s = p;

    for ( k = m1; k <= m2; k++ )
    {
      r = r + 1;
      its = 1;
      w[r-1] = rv5[k];
      x1 = rv5[k];
/*
  Look for close or coincident roots.
*/
      if ( k != m1 )
      {
        if ( eps2 <= x1 - x0 )
        {
          group = -1;
        }
        group = group + 1;
        if ( x1 <= x0 )
        {
          x1 = x0 + eps3;
        }
      }
/*
  Elimination with interchanges and initialization of vector.
*/
      v = 0.0;

      for ( i = p; i <= q; i++ )
      {
        rv6[i] = uk;

        if ( i != p )
        {
          if ( fabs ( u ) <= fabs ( e[i] ) )
          {
            xu = u / e[i];
            rv4[i] = xu;
            rv1[i-1] = e[i];
            rv2[i-1] = d[i] - x1;
            rv3[i-1] = 0.0;
            if ( i != q )
            {
              rv3[i-1] = e[i+1];
            }
            u = v - xu * rv2[i-1];
            v = - xu * rv3[i-1];
            continue;
          }

          xu = e[i] / u;
          rv4[i] = xu;
          rv1[i-1] = u;
          rv2[i-1] = v;
          rv3[i-1] = 0.0;
        }

        u = d[i] - x1 - xu * v;

        if ( i != q )
        {
          v = e[i+1];
        }
       }

      if ( u != 0.0 )
      {
        rv1[q] = u;
      }
      else
      {
        rv1[q] = eps3;
      }

      rv2[q] = 0.0;
      rv3[q] = 0.0;
/*
  Back substitution.
*/
      while ( true )
      {
        for ( i = q; p <= i; i-- )
        {
          rv6[i] = ( rv6[i] - u * rv2[i] - v * rv3[i] ) / rv1[i];
          v = u;
          u = rv6[i];
        }
/*
  Orthogonalize with respect to previous members of group.
  R and GROUP are 1-based.
*/
        for ( j = r - group; j < r - 1; j++ )
        {
          xu = 0.0;
          for ( i = p; i <= q; i++ )
          {
            xu = xu + rv6[i] * z[i+j*n];
          }
          for ( i = p; i <= q; i++ )
          {
            rv6[i] = rv6[i] - xu * z[i+j*n];
          }
        }

        norm = 0.0;
        for ( i = p; i <= q; i++ )
        {
          norm = norm + fabs ( rv6[i] );
        }

        if ( 1.0 <= norm )
        {
          break;
        }
/*
  Forward substitution.
*/
        if ( its == 5 )
        {
          ierr = 4 * n + r;
          free ( rv1 );
          free ( rv2 );
          free ( rv3 );
          free ( rv4 );
          free ( rv5 );
          free ( rv6 );
          return ierr;
        }

        if ( norm == 0.0 )
        {
          rv6[s-1] = eps4;
          s = s + 1;
          if ( q + 1 < s )
          {
            s = p + 1;
          }
        }
        else
        {
          xu = eps4 / norm;
          for ( i = p; i <= q; i++ )
          {
            rv6[i] = rv6[i] * xu;
          }
        }
/*
  Elimination operations on next vector iterate.

  If rv1[i-1] == e[i], a row interchange was performed earlier in the
  triangularization process.
*/
        for ( i = p + 1; i <= q; i++ )
        {
          u = rv6[i];

          if ( rv1[i-1] == e[i] )
          {
            u = rv6[i-1];
            rv6[i-1] = rv6[i];
          } 
          rv6[i] = u - rv4[i] * rv6[i-1];
        }
        its = its + 1;
      }
/*
  Normalize so that sum of squares is 1 and expand to full order.
*/
      u = 0.0;
      for ( i = p; i <= q; i++ )
      {
        u = pythag ( u, rv6[i] );
      }

      xu = 1.0 / u;

      for ( i = 0; i < n; i++ )
      {
        z[i+(r-1)*n] = 0.0;
      }
      for ( i = p; i <= q; i++ )
      {
        z[i+(r-1)*n] = rv6[i] * xu;
      }

      x0 = x1;
    }
  }

  free ( rv1 );
  free ( rv2 );
  free ( rv3 );
  free ( rv4 );
  free ( rv5 );
  free ( rv6 );

  return ierr;
}
