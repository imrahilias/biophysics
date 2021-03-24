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
function [ rotation, translation, s ] = rig( p, q )
  assert( nargin == 2 && size( p ) == size( q ), ...
          "need 2 identical input matrices\n" );
  assert( size( p, 2 ) == 3, "input matrix p must be nx3\n" );
  assert( size( q, 2 ) == 3, "input matrix q must be nx3\n" );
  n = size( p, 1 );
  assert( n > 2, "need at least 3 points\n" );
  
  ## find centroids & shift there
  centroid_p = mean( p, 1 );
  centroid_q = mean( q, 1 );
  p_shifted = p - repmat( centroid_p, n, 1 );
  q_shifted = q - repmat( centroid_q, n, 1 );
  
  ## covariance matrix
  h = p_shifted' * q_shifted;
  
  ## singular value decomposition
  [ u, s, v ] = svd( h );
  rotation = v * u'
  
  if ( det( rotation ) < 0 )
    printf( "warning: det(r) < 0\n" );
    if ( any( s( : ) ) < 0 )
      printf( "found reflection, correcting.\n" );
      v( :, 3 ) = -v( :, 3 );
      rotation = v * u';
    else
      printf( "error: single-value-decomposition failed! \
provided data seems is too noisy for least-squares\n." );  
    endif
  endif

  translation = centroid_q - ( rotation * centroid_p' )';
end
