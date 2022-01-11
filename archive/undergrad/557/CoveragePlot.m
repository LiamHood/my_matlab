clear ; close all ; clc ;
% ycov28 = readmatrix( 'YearCoverage28.txt' ) ;
% ycov32 = readmatrix( 'YearCoverage32.txt' ) ;
% save( 'YearCoverageConstellation.mat' ) ;
load( 'YearCoverageConstellation.mat' ) ;

figure
plot( ycov28(:,1) , ycov28(:,2) )
title( '28 Satellites' )
xlabel( 'Time Since Nominal Orbits [Days]' )
ylabel( 'Percent of Globe Uncovered' )

figure
plot( ycov32(:,1) , ycov32(:,2) )
title( '32 Satellites' )
xlabel( 'Time Since Nominal Orbits [Days]' )
ylabel( 'Percent of Globe Uncovered' )