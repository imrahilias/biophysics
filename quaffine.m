#!/bin/octave
## @ moritz siegel

clear all
close all
clc
graphics_toolkit('gnuplot')
#graphics_toolkit('qt')
#graphics_toolkit('fltk')

disp( 'init \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );

## init first point cloud.
n = 100;
theta = 100 * rand( n, 1 );
p = [ cos( theta ), sin( theta ), 1 * rand( size( theta ) ) ];

## define affine transforms & derive second point cloud.
translate = [ 0.5, 0, 0.5 ];
reflect = [ -1, 0, 0; 0, 1, 0; 0, 0, 1 ];
scale = [ 2, 0, 0; 0, 1, 0; 0, 0, 1 ];
theta = 30;
rot = [ cosd(theta), 0, -sind(theta); 0, 1, 0; sind(theta), 0, cosd(theta) ];
shear = [ 1, 0.5, 0; 0, 1, 0; 0, 0, 1 ];
affine = eye( 3 );
affine = affine * rot;
#affine = affine * reflect;
#affine = affine * scale;
affine = affine * shear;
q = ( affine * p' )';
q = q + translate;

## check determinant.
assert( det( affine ) > 0, "det( affine ) <= 0, rerun" );


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
print( fhc0, "quaffine_1_shifted_vs_origninal.png" )

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
print( fhos0, "quaffine_2_orthogonal_vs_shifted.png" );
fhoo0 = figure;
plot3( p_orthogonal(:,1), p_orthogonal(:,2), p_orthogonal(:,3), 'b' );
hold on;
axis equal;
plot3( q_orthogonal(:,1), q_orthogonal(:,2), q_orthogonal(:,3), 'r' );
print( fhoo0, "quaffine_2_orthogonal.png" );

## 3) enter 4th dimension: quaternions.
radii3_p = sum( p_orthogonal .^ 2, 2 ) / 3; # we define re( quat ) = r/3, why not unit quat?
radii3_q = sum( q_orthogonal .^ 2, 2 ) / 3;
quaternion_p = [ radii3_p, p_orthogonal ];
quaternion_q = [ radii3_q, p_orthogonal ];

## a) first-order elementary symmetric polynomials for all quaterions,
## since point clouds are orthogonal, these should be zero.
assert( sum( p_shifted ) < eps & sum( p_shifted ) < eps, "not orthogonal" );
troika_p_1 = 0;
troika_q_1 = 0;

## b) second-order elementary symmetric polynomials,
## (q1 * q2)_xyz = cross( q1_xyz, q2_xyz ) + q1_r * q2_xyz + q2_r * q1_xyz,
## (q1 * q2)_r = q1_r * q2_r - dot( q1_xyz * q2_xyz ),
## cross( a, a ) = 0, dot product is kummutative.
## https://www.sciencedirect.com/topics/computer-science/quaternion-multiplication
moment_2_p( :,2:4 ) = 2 * radii3_p .* p_orthogonal;
moment_2_q( :,2:4 ) = 2 * radii3_q .* q_orthogonal;
moment_2_p( :,1 ) = radii3_p .^ 2 - dot( p_orthogonal, p_orthogonal, 2 );
moment_2_q( :,1 ) = radii3_q .^ 2 - dot( q_orthogonal, q_orthogonal, 2 ) ;

## second-order elementary symmetric polynomials is equal to,
## minus the second-order discrete moment for all quaterions.
troika_p_2 = - sum( moment_2_p );
troika_q_2 = - sum( moment_2_q );

## c) third-order elementary symmetric polynomials,
## which is equal to the third-order discrete moment for all quaterions.
moment_3_p( :,2:4 ) = moment_2_p( :,1 ) .* p_orthogonal + radii3_p .* moment_2_p( :,2:4 );
moment_3_q( :,2:4 ) = moment_2_q( :,1 ) .* q_orthogonal + radii3_q .* moment_2_q( :,2:4 );
moment_3_p( :,1 ) = moment_2_p( :,1 ) .* radii3_p - dot( moment_2_p( :,2:4 ), p_orthogonal, 2 );
moment_3_q( :,1 ) = moment_2_q( :,1 ) .* radii3_q - dot( moment_2_q( :,2:4 ), q_orthogonal, 2 );
troika_p_3 = moment_3_p;
troika_q_3 = moment_3_q;

## horns algorithm s_ab = sum_k( p_a_k * p_b_k ), here k=(1),2,3; since troika_q_1 = 0
s_xx = troika_p_2( 1 ) * troika_q_2( 1 ) + troika_p_3( 1 ) * troika_q_3( 1 );
s_xy = troika_p_2( 1 ) * troika_q_2( 2 ) + troika_p_3( 1 ) * troika_q_3( 2 );
s_xz = troika_p_2( 1 ) * troika_q_2( 3 ) + troika_p_3( 1 ) * troika_q_3( 3 );
s_yx = troika_p_2( 2 ) * troika_q_2( 1 ) + troika_p_3( 2 ) * troika_q_3( 1 );
s_yy = troika_p_2( 2 ) * troika_q_2( 2 ) + troika_p_3( 2 ) * troika_q_3( 2 );
s_yz = troika_p_2( 2 ) * troika_q_2( 3 ) + troika_p_3( 2 ) * troika_q_3( 3 );
s_zx = troika_p_2( 3 ) * troika_q_2( 1 ) + troika_p_3( 3 ) * troika_q_3( 1 );
s_zy = troika_p_2( 3 ) * troika_q_2( 2 ) + troika_p_3( 3 ) * troika_q_3( 2 );
s_zz = troika_p_2( 3 ) * troika_q_2( 3 ) + troika_p_3( 3 ) * troika_q_3( 3 );

## matrix of sums of products (symetric).
## note that the trace, sum of diagonal elements of the matrix n, is zero.
## this takes care of the 10th degree of freedom.
big_n( 1, 1 ) = s_xx + s_yy + s_zz;
big_n( 1, 2 ) = s_yz - s_zy;
big_n( 1, 3 ) = s_zx - s_xz;
big_n( 1, 4 ) = s_xy - s_yx;
big_n( 2, 1 ) = big_n( 1, 2 );
big_n( 2, 2 ) = s_xx - s_yy - s_zz;
big_n( 2, 3 ) = s_xy + s_yx;
big_n( 2, 4 ) = s_zx + s_xz;
big_n( 3, 1 ) = big_n( 1, 3 );
big_n( 3, 2 ) = big_n( 2, 3 );
big_n( 3, 3 ) = - s_xx + s_yy - s_zz;
big_n( 3, 4 ) = s_yz + s_zy;
big_n( 4, 1 ) = big_n( 1, 4 );
big_n( 4, 2 ) = big_n( 2, 4 );
big_n( 4, 3 ) = big_n( 3, 4 );
big_n( 4, 4 ) = - s_xx - s_yy + s_zz;

## the unit quaternion (rotation) that maximizes,is the eigenvector
## corresponding to the most positive eigenvalue of the matrix N.
## solve eigenvalue problem return eigenvector correspodning to the
## largest (algebraic) eigenvalue lambda.
[ v, lambda ] = eigs ( big_n, 1, "la" ); # v not unit quaternion!

## construct a rotation matrix out of the quaternion v.
q0 = v(1);
q1 = v(2);
q2 = v(3);
q3 = v(4);
rotation( 1,1 ) = 2 * (q0 * q0 + q1 * q1) - 1;
rotation( 1,2 ) = 2 * (q1 * q2 - q0 * q3);
rotation( 1,3 ) = 2 * (q1 * q3 + q0 * q2);
rotation( 2,1 ) = 2 * (q1 * q2 + q0 * q3);
rotation( 2,2 ) = 2 * (q0 * q0 + q2 * q2) - 1;
rotation( 2,3 ) = 2 * (q2 * q3 - q0 * q1);
rotation( 3, 1 ) = 2 * (q1 * q3 - q0 * q2);
rotation( 3, 2 ) = 2 * (q2 * q3 + q0 * q1);
rotation( 3, 3 ) = 2 * (q0 * q0 + q3 * q3) - 1;

## transform the rotation in the non-orthogonal affine transform.
affine_rotation = ( s_sqrt_q * rotation ) * s_inv_sqrt_p;

## rotate orthogonal point clouds.
p_rotated = ( rotation * p_orthogonal' )';

## plot & compare both orthonormalised point clouds with shifted.
fhor0 = figure;
plot3( p_orthogonal(:,1), p_orthogonal(:,2), p_orthogonal(:,3), '.b' );
hold on;
axis equal;
plot3( q_orthogonal(:,1), q_orthogonal(:,2), q_orthogonal(:,3), 'r' );
plot3( p_rotated(:,1), p_rotated(:,2), p_rotated(:,3), 'b' );
print( fhor0, "quaffine_2_rotated_vs_orthogonal.png" );

## it is, unfortunately not working!
