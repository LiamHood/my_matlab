function dstate_dt = TwoBodyMotion( t , state , mu ) 
% Finds change in state with respect to time. Input time, t, in seconds and
% state as position vector followed by velocity vector as well as mu

rad = norm( [ state(1) state(2) state(3) ] ) ; %radius

dx = state(4) ; % velocity in x
dy = state(5) ; % velocity in y
dz = state(6) ; % velocity in z
ddx = -mu*state(1)/rad^3 ; % acceleration in x
ddy = -mu*state(2)/rad^3 ; % acceleration in y
ddz = -mu*state(3)/rad^3 ; % acceleration in z

dstate_dt = [ dx ; dy ; dz ; ddx ; ddy ; ddz ] ;

end