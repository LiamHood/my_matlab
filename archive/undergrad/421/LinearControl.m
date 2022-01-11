function [ dstate , t ] = LinearControl( t , state , I , kp , kd , ceta , ceps )
eps(1:3,1) = state( 1:3 ) ;
eta = state( 4 ) ;
quat = [ eta ; eps ]' ;
cquat = [ ceta ; ceps ]' ;
qstar = quatconj( cquat ) ;
qerr = quatmultiply( qstar , quat ) ;
etae = qerr( 1 ) ;
epse(1:3,1) = qerr( 2:4 ) ;
w(1:3,1) = state( 5:7 ) ;

deps = w./2 ;
deta = 0 ;
dw = inv( I )*( -kp*epse - kd*w ) ;

dstate = [ deps ; deta ; dw ] ;

end