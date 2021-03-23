#!/bin/octave
## @ moritz siegel

clear all
close all
clc
graphics_toolkit('gnuplot')
#graphics_toolkit('qt')
#graphics_toolkit('fltk')

disp( 'init \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );

## init first point cloud
n = 100;
theta = 100 * rand( n, 1 );
p = [ cos( theta ), sin( theta ), 1 * rand( size( theta ) ) ];

## define affine transforms & derive second point cloud
translate = [ 0.5, 0, 0.5 ];
reflect = [ -1, 0, 0; 0, 1, 0; 0, 0, 1 ];
scale = [ 2, 0, 0; 0, 1, 0; 0, 0, 1 ];
theta = 30;
rot = [ cosd(theta), 0, -sind(theta); 0, 1, 0; sind(theta), 0, cosd(theta) ];
shear = [ 1, 0.5, 0; 0, 1, 0; 0, 0, 1 ];
affine = eye( 3 );
affine = affine * rot;
#affine = affine * reflect;
affine = affine * scale;
affine = affine * shear;
q = ( affine * p' )';
q = q + translate;

## check determinant
assert( det( affine ) > 0, "det( affine ) <= 0, rerun" );


disp( 'analyse \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );

## 1) get rid of translation
centroid_p = mean( p );
centroid_q = mean( q );
p_shifted = p - repmat( centroid_p, n, 1 );
q_shifted = q - repmat( centroid_q, n, 1 );

## plot & compare both shifted point clouds with originals
fhc0 = figure;
plot3( p(:,1), p(:,2), p(:,3), '.b' );
hold on
plot3( q(:,1), q(:,2), q(:,3), '.r' );
axis equal
plot3( p_shifted(:,1), p_shifted(:,2), p_shifted(:,3), 'b' );
plot3( q_shifted(:,1), q_shifted(:,2), q_shifted(:,3), 'r' );
print( fhc0, "quaffine_1_shifted_vs_origninal.png" )

## check centroids
eps = 1e-8;
assert( mean( p_shifted ) < eps & mean( q_shifted ) < eps, "not centered" );

## 2) orthogonal reduction

## covariance matrices
s_p = p_shifted' * p_shifted;
s_q = q_shifted' * q_shifted;
#assert( [1,0,0] * s_p * [1,0,0]' > 0 & [1,0,0] * s_q * [1,0,0]' > 0, "not positive-definite"); # chol() checks that

## inverse square-roots of the covariance matrices
l_p = chol( s_p, "lower" ); # choleski decomposition: A -> LL*
l_q = chol( s_q, "lower" );
ss_p = inv( l_p ); # inverse: A -> A^-1
ss_q = inv( l_q );
p_orthogonal = ( ss_p * p_shifted' )';
q_orthogonal = ( ss_q * q_shifted' )';

## plot & compare both orthonormalised point clouds with shifted
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

