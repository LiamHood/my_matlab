function [ dstate , t ] = NonLinearControlne( t , state , I , kp , kd )
eps(1:3,1) = state( 1:3 ) ;
eta = state( 4 ) ;
w(1:3,1) = state( 5:7 ) ;
wcross = [ 0 -w(3) w(2) ; w(3) 0 -w(1) ; -w(2) w(1) 0 ] ;
epscross = [ 0 -eps(3) eps(2) ; eps(3) 0 -eps(1) ; -eps(2) eps(1) 0 ] ;

T = -kp*eps - kd*w ;

deps = .5*( eta*eye( 3 ) + epscross )*w ;
deta = -.5*eps'*w ;
dw = inv( I )*( -wcross*I*w + T ) ;

dstate = [ deps ; deta ; dw ] ;

end