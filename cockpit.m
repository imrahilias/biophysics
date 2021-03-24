#!/bin/octave
## images fluorescence microscopy / 3d position / offsets.
## splice red/green/blank fluorescence microscopy images already located
## into red/green matrices for further analysis, find rigid transform,
## do statistics on all files.
## @ moritz siegel

## init
close all
clear all
clc

## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## soft settings \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
global hd = "~/biophysics" # parent directory
global wd = "data/210318_beads_2_colors_in_focus"
global fn = "position"
global nn = ""; # additional luminosity suffix  ( default: "" )
method = "affine"; # /char, "rigid","affine", method to reconstruct
global minpts = 10; # minimum number of points needed in its neighbourhood to consider it as a valid data(not noise). ( default: 3 )
global dist = 300; #/nm, distance on which neighbourhood is calculated.( default: 200 )
nmax = 1; #/num, maximum number of files analysed.( default: 10 )
global verbose = true; # /bool, plot for error checking? ( default: true )
exclude = [ 6, 7, 11 ]
channel_order0 = [ 0, 1, 2 ];
channel_order = repmat( channel_order0, nmax, 1 );
channel_order( 3, : ) = [ 2, 0, 1 ];
channel_order( 4, : ) = [ 1, 2, 0 ];
channel_order( 6, : ) = [ 2, 0, 1 ];
channel_order( 7, : ) = [ 1, 2, 0 ];
channel_order( 8, : ) = [ 2, 0, 1 ];
channel_order( 12, : ) = [ 2, 0, 1 ];
## end of settings: hands off /////////////////////////////////////////////////
## ///////////////////////////////////////////////////////////////////////////


disp( 'warmup \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );

## lets move it move it
addpath( hd ) # function & stuff needed
stamp = strftime("%Y_%m_%d_%H%M%S", localtime (time ())); # create timestamp for saving files
global nwd = sprintf( "%s/%s/all_%s%s", hd, wd, stamp, nn ); # concat working dir
mkdir( nwd )
chdir( nwd )
global rf = fopen ( "kockpit.log", "w" );


disp( 'cycle \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );
rotations = nan( nmax, 3, 3 );
translations = nan( nmax, 3 );
for n = 1 : nmax
    if ( any( exclude == n ) )
        disp( sprintf( "    warning: skipping file positions%d.csv", n ) );
        continue
    endif
    disp( sprintf( "    importing positions%d.csv", n ) );
    disp( '    load \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );
    
    ## "id","frame","x [nm]","y [nm]","sigma [nm]","intensity [photon]","offset [photon]","bkgstd [photon]","uncertainty [nm]"
    try
        pos = csvread( sprintf( "%s/%s/%s%d.csv", hd, wd, fn, n ) );
    catch
        continue
    end_try_catch
    
    
    disp( '    analyse \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );
    
    [ centroid_red, centroid_blue_transformed, centroid_blue, red, blue, rotation, translation ] = ampel( pos, channel_order( n, : ), n, method );
    
    ## save optimal rigid transform for comparison
    rotations( n, :, : ) = rotation;
    translations( n, : ) = translation;
    
    
    disp( 'plot \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );
    
    if ( verbose  == true )
        if ( n == 1 )
            fh0 = figure;
        else
            clear fh0
        endif     
        fh1 = scatter( blue( :, 3 ), blue( :, 4 ), '.b');
        hold on
        fh2 = scatter( red( :, 3 ), red( :, 4 ), '.r');
        fh3 = plot3( centroid_blue(:,1),...
        centroid_blue(:,2),...
        centroid_blue(:,3), 'bx' );
        fh4 = plot3( centroid_red(:,1),...
        centroid_red(:,2),...
        centroid_red(:,3), 'rx' );
        fh5 = plot3( centroid_blue_transformed(:,1),...
        centroid_blue_transformed(:,2), ...
        centroid_blue_transformed(:,3), 'b+' );
    endif
    
endfor

# remove empty/corrupt lines for statistics
n = 1
while ( n <= nmax )
    if ( any( isnan( rotations( n, :, : ) ) ) )
        rotations( n, :, : ) = []
        translations( n, : ) = []
        nmax = nmax - 1;
    else
        n = n + 1;
    endif
endwhile

save rotations.mat rotations
save translations.mat translations

mean( rotations )
std( rotations )
median( rotations )
figure;
plot3( rotations(:,1), rotations(:,2), rotations(:,3) )

mean( translations )
std( translations )
median( translations )
plot3( translations(:,1), translations(:,2), translations(:,3) )


disp( 'done \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );
