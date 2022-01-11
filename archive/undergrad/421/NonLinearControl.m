function [ dstate , t ] = NonLinearControl( t , state , I , kp , kd , ceta , ceps )
eps(1:3,1) = state( 1:3 ) ;
eta = state( 4 ) ;
quat = [ eta ; eps ]' ;
cquat = [ ceta ; ceps ]' ;
qstar = quatconj( cquat ) ;
qerr = quatmultiply( qstar , quat ) ;
etae = qerr( 1 ) ;
epse(1:3,1) = qerr( 2:4 ) ;
w(1:3,1) = state( 5:7 ) ;

w(1:3,1) = state( 5:7 ) ;
wcross = [ 0 -w(3) w(2) ; w(3) 0 -w(1) ; -w(2) w(1) 0 ] ;
epscrosse = [ 0 -epse(3) epse(2) ; epse(3) 0 -epse(1) ; -epse(2) epse(1) 0 ] ;

T = -kp*epse - kd*w ;

deps = .5*( etae*eye( 3 ) + epscrosse )*w ;
deta = -.5*epse'*w ;
dw = inv( I )*( -wcross*I*w + T ) ;

dstate = [ deps ; deta ; dw ] ;

end