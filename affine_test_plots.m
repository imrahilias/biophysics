#!/bin/octave
## @ moritz siegel

clear all
close all
clc
#graphics_toolkit('gnuplot')
#graphics_toolkit('qt')
#graphics_toolkit('fltk')

disp( 'init \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );

## init first point cloud.
n = 100;
theta = 100 * rand( n, 1 );
p = [ cos( theta ), sin( theta ), 1 * rand( size( theta ) ) ];

## define affine transforms & derive second point cloud.
translate = [ 0, 0.3, 0.7 ];
reflect = [ -1, 0, 0; 0, 1, 0; 0, 0, 1 ];
scale = [ 2, 0, 0; 0, 1, 0; 0, 0, 1 ];
theta = 30;
rot = [ cosd(theta), 0, -sind(theta); 0, 1, 0; sind(theta), 0, cosd(theta) ];
shear = [ 1, 0.5, 0; 0, 1, 0; 0, 0, 1 ];
affine_transform = eye( 3 );
affine_transform = affine_transform * rot;
affine_transform = affine_transform * reflect;
affine_transform = affine_transform * scale;
affine_transform = affine_transform * shear;
q = ( affine_transform * p' )';
q = q + translate;

## check determinant.
#assert( det( affine_transform ) > 0, "det( affine_transform ) <= 0, rerun" );


disp( 'analyse \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );

## 1) get rid of translation.
centroid_p = mean( p );
centroid_q = mean( q );
p_shifted = p - repmat( centroid_p, n, 1 );
q_shifted = q - repmat( centroid_q, n, 1 );

## plot & compare both shifted point clouds with originals.
fhc0 = figure;
plot3( p(:,1), p(:,2), p(:,3), '.b' );
hold on
plot3( q(:,1), q(:,2), q(:,3), '.r' );
axis equal
plot3( p_shifted(:,1), p_shifted(:,2), p_shifted(:,3), 'b' );
plot3( q_shifted(:,1), q_shifted(:,2), q_shifted(:,3), 'r' );
print( fhc0, "affine_1_shifted_vs_origninal.png" )

## check centroids.
eps = 1e-8;
assert( mean( p_shifted ) < eps & mean( q_shifted ) < eps, "not centered" );

## 2) orthogonal reduction

## covariance matrices.
s_p = p_shifted' * p_shifted;
s_q = q_shifted' * q_shifted;
#assert( [1,0,0] * s_p * [1,0,0]' > 0 & [1,0,0] * s_q * [1,0,0]' > 0, "not positive-definite"); # chol() checks that

## inverse square-roots of the covariance matrices.
s_sqrt_p = chol( s_p, "lower" ); # choleski decomposition: A -> LL*
s_sqrt_q = chol( s_q, "lower" );
s_inv_sqrt_p = inv( s_sqrt_p ); # inverse: A -> A^-1
s_inv_sqrt_q = inv( s_sqrt_q );
p_orthogonal = ( s_inv_sqrt_p * p_shifted' )';
q_orthogonal = ( s_inv_sqrt_q * q_shifted' )';

## p and q are now related by a rotational matrix r = s_inv_sqrt_q affine s_sqrt_p (no inverse!)

## plot & compare both orthonormalised point clouds with shifted.
fhos0 = figure;
plot3( p_shifted(:,1), p_shifted(:,2), p_shifted(:,3), '.b' );
hold on;
axis equal;
plot3( q_shifted(:,1), q_shifted(:,2), q_shifted(:,3), '.r' );
plot3( p_orthogonal(:,1), p_orthogonal(:,2), p_orthogonal(:,3), 'b' );
plot3( q_orthogonal(:,1), q_orthogonal(:,2), q_orthogonal(:,3), 'r' );
print( fhos0, "affine_2_orthogonal_vs_shifted.png" );
fhoo0 = figure;
plot3( p_orthogonal(:,1), p_orthogonal(:,2), p_orthogonal(:,3), 'b' );
hold on;
axis equal;
plot3( q_orthogonal(:,1), q_orthogonal(:,2), q_orthogonal(:,3), 'r' );
print( fhoo0, "affine_3_orthogonal.png" );

## covariance matrix
h = p_orthogonal' * q_orthogonal;

## singular value decomposition
[ u, s, v ] = svd( h );
rotation = v * u';

if ( det( rotation ) < 0 )
  printf( "warning: det(rotation) < 0" );
  if ( any( s(:) ) < 0 )
    printf( "found reflection, correcting." );
    v(:,3) = -v(:,3);
    rotation = v*u';
  else
    printf( "error: single-value-decomposition failed! \
provided data seems is too noisy for least-squares." );  
  endif
endif

## rotate orthogonal point clouds.
p_rotated = ( rotation * p_orthogonal' )';

## plot & compare both orthonormalised point clouds with shifted.
fhor0 = figure;
plot3( p_orthogonal(:,1), p_orthogonal(:,2), p_orthogonal(:,3), '.b' );
hold on;
axis equal;
plot3( q_orthogonal(:,1), q_orthogonal(:,2), q_orthogonal(:,3), '.r' );
plot3( p_rotated(:,1), p_rotated(:,2), p_rotated(:,3), 'b' );
print( fhor0, "affine_4_rotated_vs_orthogonal.png" );

## transform the rotation in the non-orthogonal affine transform.
affine_rotation = ( s_sqrt_q * rotation ) * s_inv_sqrt_p;

## rotate orthogonal point clouds.
p_shifted_rotated = ( affine_rotation * p_shifted' )';

## plot & compare both orthonormalised point clouds with shifted.
fhor0 = figure;
plot3( p_shifted(:,1), p_shifted(:,2), p_shifted(:,3), '.b' );
hold on;
axis equal;
plot3( q_shifted(:,1), q_shifted(:,2), q_shifted(:,3), '.r' );
plot3( p_shifted_rotated(:,1), p_shifted_rotated(:,2), p_shifted_rotated(:,3), 'b' );
print( fhor0, "affine_5_affine_rotated_vs_shifted.png" );

## compute translation.
translation = centroid_q - ( affine_rotation * centroid_p' )';

disp( 'simulated translation vector:' )
disp( translate )
disp( 'recovered translation vector:' )
disp( translation )
disp( 'deviation translation vector:' )
disp( translate - translation )
disp( 'simulated affine transform:' )
disp( affine_transform )
disp( 'recovered affine transform:' )
disp( affine_rotation )
disp( 'deviation affine transform:' )
disp( affine_transform - affine_rotation )
assert( all( all( ( affine_transform - affine_rotation ) < eps ) ), "failed to recover affine transform" ); 
