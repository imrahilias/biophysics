#!/bin/octave
## @ moritz siegel
function[ centroid_red, centroid_blue_transformed, centroid_blue, red, blue, r, t ] = ampel( pos, channel_order, n );
## images fluorescence microscopy / 3d position / offsets.
## splice red/green/blank fluorescence microscopy images already located
## into red/green matrices for further analysis, find rigid transform.

global nwd
global rf
global minpts
global eps

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
[ assignments_blue, c_blue ] = dbscan( blue( :, 3:4 ), minpts, eps );
centroid_blue = zeros( c_blue, 3 ); # x y z components
for k = 1 : c_blue
    idx = find( assignments_blue == k );
    centroid_blue( k, 1:2 ) = mean( blue( idx, 3:4 ), 1 );
endfor
[ assignments_red, c_red ] = dbscan( red( :, 3:4 ), minpts, eps );
centroid_red = zeros( c_red, 3 );
for k = 1 : c_red
    idx = find( assignments_red == k );
    centroid_red( k, 1:2 ) = mean( red( idx, 3:4 ), 1 );
endfor

## test affine transforms
##theta = 30;
##rot = [ cosd(theta), 0, -sind(theta); 0, 1, 0; sind(theta), 0, cosd(theta) ];
##shear = [ 1, 0.5, 0; 0, 1, 0; 0, 0, 1 ];
##scale = [ 2, 0, 0; 0, 1, 0; 0, 0, 1 ];
###affine = rot * shear * scale;
##affine = rot;
##centroid_red = ( affine * centroid_red' )';

## find optimal rigid transform
[ r, t ] = rig( centroid_blue', centroid_red' );

## transform blue set to red set
centroid_blue_transformed = ( r * centroid_blue' + t )';

## save matrices
swd = sprintf( "%s/position%d", nwd, n ); # concat working dir
mkdir( swd )
chdir( swd )
save red.mat red
save blue.mat blue
save empty.mat empty
save centroid_red.mat centroid_red
save centroid_blue.mat centroid_blue
save r.mat r
save t.mat t
save centroid_blue_transformed.mat centroid_blue_transformed

endfunction
