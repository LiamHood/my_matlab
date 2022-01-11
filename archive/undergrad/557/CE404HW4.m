clear ; close all ; clc ;
pointx = 1.3229e+006 ;
pointy = -8.9239e+005 ;
pointx2 = -1.0502e+006 ;

center = ( pointx + pointx2 )/2 ;
taumax = sqrt( ( pointx - center )^2 + pointy^2 ) ;

x = linspace( center - taumax , center + taumax , 1e3 ) ;
y1 = sqrt( taumax^2 - (x-center).^2 ) ;
y2 = -sqrt( taumax^2 - (x-center).^2 ) ;



figure
hold on
axis square
plot( pointx , pointy , '*r' )
plot( x , y1 , 'k' )
plot( x , y2 , 'b' )
plot( pointx2 , -pointy , '*r' )
legend( 'Stress State at Point' , 'Mohr''s Circle' , 'Location' , 'northeast' )
xlabel( 'Normal Stress [Pa]' )
ylabel( 'Shear Stress [Pa]' )
title( 'Stress at Point of Maximum Stress Normal to Weld Line' )