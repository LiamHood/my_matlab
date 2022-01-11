function [ dstate , t ] = TwoBody( t , state , mu )
    r = state(1:3) ;
    v = state(4:6) ;
    rad = norm( [ r(1) r(2) r(3) ] ) ; %radius
    dstate = zeros(6,1) ;
    dstate(1) = v(1) ; % velocity in x
    dstate(2) = v(2) ; % velocity in y
    dstate(3) = v(3) ; % velocity in z
    dstate(4) = -mu*r(1)/rad^3 ; % acceleration in x
    dstate(5) = -mu*r(2)/rad^3 ; % acceleration in y
    dstate(6) = -mu*r(3)/rad^3 ; % acceleration in z
end