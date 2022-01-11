clear ; close all ; clc ;

% results = xlsread( 'Tuetilever.csv' ) ;
% save( 'Cantilever.mat' ) ;

load( 'Cantilever.mat' ) ;
P = 4 ;
E = 73.1e9 ;
L = .3 ;
density = 2.7e-3*1e6 ;
A = 60.57e-3*1.72e-3 ;

strain = results( : , 4 ) ;
x = results( : , 5 )*1e-3 ;
weight = density*A*9.8 ;
c = results( : , 9 )*1e-3 ;
I = 25.44e-12 ;
M = P.*(L-x) ;
Mw = weight.*(L-x).^2/2 ;
DeflectionTHEO = (P/(6*E*I))*(x.^3-3*L*x.^2) ;
Deflection = results( : , 6 ) ;
dmean = mean( Deflection ) ;
dstd = std( Deflection ) ;
dTheo = P*L^3/(3*E*I)*1e3 ;



figure
hold on
plot( x*1e2 , M , '*b'  )
% plot( x*1e2 , Mw+M , '*r' )
% legend( 'With only applied load' , 'With weight of beam' )
ylabel( 'Bending Moment [Nm]' )
xlabel( 'Distance from Support [cm]' )
hold off




