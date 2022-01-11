function [ dstate ] = Qrates( t , state , w0 )
% Find change in quaternions using current quaternion position and angular
% velocity with no torque
eps = state( 1:3 ) ;
eta = state( 4 ) ;
C = quat2rotm( eta , eps ) ;
w = C*w0 ;
epscross = crossmatrix( eps ) ;
epsdot = .5 * ( eta*eye(3) + epscross ) * w ;
etadot = .5 * eps' * w ;
dstate = [ epsdot' , etadot ]' ;
    function [ across ] = crossmatrix( a )
        across = [ 0 -a(3) a(2) ; ...
                   a(3) 0 -a(1) ; ...
                  -a(2) a(3) 0 ] ;
    end
end
