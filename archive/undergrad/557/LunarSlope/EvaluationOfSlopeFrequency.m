clear ; close all ; clc ;
load( 'slopeFrequency.mat' )
% load( 'MaxSlopeFrequency.mat' )
% load( 'AvgAbsSlopeFrequency.mat' )
load( 'AvgSlopeFrequency.mat' )
NslopeFreq = slopeFreq./( sum( slopeFreq ) ) ;
for ii = 1:length( slopeFreq )
    percentSlopes(ii) = sum( NslopeFreq( 1:ii ) ) ;
end

figure
plot( posSlopes , NslopeFreq )
xlabel( 'Necessary Elevation Angle' )
ylabel( 'Normalized amount of points' )
title( 'Elevation angle required to clear terrain' )

I = find( posSlopes == 16 ) ;
figure
hold on
plot( posSlopes , percentSlopes*100 )
plot( posSlopes(I) , percentSlopes(I)*100 , '.r' )
xlabel( 'Minimum Elevation Angle' )
ylabel( 'Percent of Points Covered' )
