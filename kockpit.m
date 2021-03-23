#!/bin/octave
## images fluorescence microscopy / 3d position / offsets.
## splice red/green/blank fluorescence microscopy images already located
## into red/green matrices for further analysis, find rigid transform,
## do statistics on all files.
## @ moritz siegel

## init
pkg load matgeom
close all
clear all
clc

## \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
## soft settings \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
global hd = "~/biophysics" # parent directory
global wd = "data/210318_beads_2_colors_in_focus"
global fn = "position"
global nn = ""; # additional luminosity suffix  ( default: "" )
global minpts = 10; # minimum number of points needed in its neighbourhood to consider it as a valid data(not noise). ( default: 3 )
global eps = 300; #/nm, distance on which neighbourhood is calculated.( default: 200 )
nmax = 15; #/num, maximum number of files analysed.( default: 10 )
global verbose = true; # /bool, plot for error checking? ( default: true )
exclude = [ 6, 7, 11 ]
channel_order0 = [ 0, 1, 2 ];
channel_order = repmat( channel_order0, nmax, 1 )
channel_order( 3, : ) = [ 2, 0, 1 ]
channel_order( 4, : ) = [ 1, 2, 0 ]
channel_order( 6, : ) = [ 2, 0, 1 ]
channel_order( 7, : ) = [ 1, 2, 0 ]
channel_order( 8, : ) = [ 2, 0, 1 ]
channel_order( 12, : ) = [ 2, 0, 1 ]
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
rs = nan( nmax, 3, 3 );
ts = nan( nmax, 3 );
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
    
    [ centroid_red, centroid_blue_transformed, centroid_blue, red, blue, r, t ] = ampel( pos, channel_order( n, : ), n );
    
    ## save optimal rigid transform for comparison
    rs( n, :, : ) = r;
    ts( n, : ) = t;
    
    
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
    if ( any( isnan( rs( n, :, : ) ) ) )
        rs( n, :, : ) = []
        ts( n, : ) = []
        nmax = nmax - 1;
    else
        n = n + 1;
    endif
endwhile

save rs.mat rs
save ts.mat ts

mean( rs )
std( rs )
median( rs )
figure;
plot3( rs(:,1), rs(:,2), rs(:,3) )

mean( ts )
std( ts )
median( ts )
plot3( ts(:,1), ts(:,2), ts(:,3) )


disp( 'done \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ ' );
