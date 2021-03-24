#!/bin/octave
## @ moritz siegel
function[ centroid_red, centroid_blue_transformed, centroid_blue, red, blue, rotation, translation ] = ampel( pos, channel_order, n, method );
## images fluorescence microscopy / 3d position / offsets.
## splice red/green/blank fluorescence microscopy images already located
## into red/green matrices for further analysis, find rigid transform.

global nwd
global rf
global minpts
global dist

## sort red/blue/empty of many frames 
red = blue = empty = zeros( size( pos ));
nb = nr = ne = 1;
for k = 1 : size( pos, 1 )
    
    ## exclude blank lines (text)
    if ( all( pos( k, : ) == 0 ) )
        disp( sprintf( 'warning: line %d is empty; skipping.', k ));
        fprintf( rf, 'warning: empty line; skipping.\n' );
        continue
    endif
    
    ## red pill / blue pill?
    channel = mod( pos( k, 2 ), 3 );
    switch ( channel )
        case channel_order(1)
            blue( nb, : ) = pos( k, : );
            nb = nb + 1;
        case channel_order(2)
            empty( ne, : ) = pos( k, : );
            ne = ne + 1;
        case channel_order(3)
            red( nr, : ) = pos( k, : );
            nr = nr + 1;
    endswitch
endfor
disp( sprintf( 'sorted localisations:\n %d red\n %d blue\n %d empty', nr, nb, ne ) );
fprintf( rf, 'sorted localisations:\n %d red\n %d blue\n %d empty\n', nr, nb, ne );

## analyse clusters & average centroids
[ assignments_blue, c_blue ] = dbscan( blue( 1:300, 3:4 ), minpts, dist );
for k = 1 : c_blue
    idx = find( assignments_blue == k );
    centroid_blue( k, 1:2 ) = mean( blue( idx, 3:4 ), 1 );
endfor
[ assignments_red, c_red ] = dbscan( red( 1:300, 3:4 ), minpts, dist );
for k = 1 : c_red
    idx = find( assignments_red == k );
    centroid_red( k, 1:2 ) = mean( red( idx, 3:4 ), 1 );
endfor

## need to ditch clusters if the sets are not equally sized
assert( c_red == c_blue, "clusters are not equally sized" );

## recover transform
switch( method )
  case "rigid"
    centroid_red = [ centroid_red, zeros( c_red, 1 ) ] # we only have 2d data currently
    centroid_blue = [ centroid_blue, zeros( c_blue, 1 ) ]
    [ rotation, translation ] = rig( centroid_blue, centroid_red )
case "affine"
    z = 1 + 1e-3 * rand( c_red, 1 );
    centroid_red = [ centroid_red, z ]; # we only have 2d data currently
    centroid_blue = [ centroid_blue, z ];
    [ rotation, translation ] = affine( centroid_blue, centroid_red );
endswitch

## transform blue set to red set
centroid_blue_transformed = ( rotation * centroid_blue' )' + translation;

## save matrices
swd = sprintf( "%s/position%d", nwd, n ); # concat working dir
mkdir( swd )
chdir( swd )
save red.mat red
save blue.mat blue
save empty.mat empty
save centroid_red.mat centroid_red
save centroid_blue.mat centroid_blue
save rotation.mat rotation
save translation.mat translation
save centroid_blue_transformed.mat centroid_blue_transformed

endfunction
