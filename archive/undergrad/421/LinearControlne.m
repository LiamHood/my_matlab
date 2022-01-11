function [ dstate , t ] = LinearControlne( t , state , I , kp , kd )
eps(1:3,1) = state( 1:3 ) ;
eta = state( 4 ) ;
w(1:3,1) = state( 5:7 ) ;

deps = w./2 ;
deta = 0 ;
dw = inv( I )*( -kp*eps - kd*w ) ;

dstate = [ deps ; deta ; dw ] ;

end