#!/bin/octave
## this function rig(p,q) finds the optimal rigid transform in 3-dimensional
## euclidian space, using least squares and
## single-value-decomposition. input a 3xn matrix (set of n 3d
## positions), and the rotation matrtix r, it returns the translation
## vector t.
## K. S. Arun, T. S. Huang and S. D. Blostein, "Least-Squares Fitting of
## Two 3-D Point Sets," in IEEE Transactions on Pattern Analysis and
## Machine Intelligence, vol. PAMI-9, no. 5, pp. 698-700, Sept. 1987,
## doi: 10.1109/TPAMI.1987.4767965.
## moritz siegel @ 210322
function [ r, t, s ] = rig( p, q )
  assert( nargin == 2 && size( p ) == size( q ), ...
          'need 2 identical input matrices' );
  assert( size( p, 1 ) == 3, 'input matrix p must be 3xn' );
  assert( size( q, 1 ) == 3, 'input matrix q must be 3xn' );
  n = size( p, 2 );
  assert( n > 2, 'need at least 2 points' );
  
  ## find centroids & shift there
  centroid_p = mean( p, 2 );
  centroid_q = mean( q, 2 );
  p_shifted = p - repmat( centroid_p, 1, n );
  q_shifted = q - repmat( centroid_q, 1, n );
  
  ## covariance matrix
  h = p_shifted * q_shifted';
  
  ## singular value decomposition
  [ u, s, v ] = svd( h );
  r = v * u';
  
  if ( det( r ) < 0 )
    printf( "warning: det(r) < 0" );
    if ( any( s(:) ) < 0 )
      printf( "found reflection, correcting." );
      v(:,3) = -v(:,3);
      r = v*u';
    else
      printf( "error: single-value-decomposition failed! \
provided data seems is too noisy for least-squares." );  
      exit(1)
    endif
  endif

  t = centroid_q - r * centroid_p;
end
