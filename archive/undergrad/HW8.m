clear ; close all ; clc ;
I = [ 1200 0 0 ; 0 2000 0 ; 0 0 2800 ] ;
Td = [ 1 1 1 ]' ;
q0 = [ .5 -.5 -.5 .5 ]' ;
eps0 = q0( 1:3 ) ;
w0 = [ 0 0 0 ] ;
ts = 30 ;
zeta = .65 ; % Damping coefficient
wn = log( 0.02*sqrt( 1 - zeta^2 ) )/( -zeta*ts ) ;
c = 2*zeta*wn ;
wd = wn*sqrt( 1 - zeta^2 ) ;
kp = 2.*I.*wn^2 ;
kd = I.*2*zeta*wn ;
t = linspace( 0 , 1e2 , 1e3 ) ;
y = 1 - exp( -zeta*wn.*t ).*( cos(wd.*t) + (zeta*wn/wd).*sin(wd.*t) ) ;
figure 
plot( t , y ) 
