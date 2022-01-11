clear ; close all ; clc ;

k = 1.4 ;
R = 1717 ;
M = linspace( .1 , 2 , 1e2 ) ;
MFP = M.*sqrt( k/R ).*( 1 + ( ( k-1 )/2 ).*M.^2 ) ;

figure
plot( M , MFP ) 
xlabel( 'Mach' )
ylabel( 'MFP [ sqrt(R)*s/ft ]' )