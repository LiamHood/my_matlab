    function dstate_dt = NonImpulsive( t , state , mu , thrust , isp , g0 ) 
    % Finds change in state with respect to time. Input time, t, in seconds and
    % state as position vector followed by velocity vector as well as mu

    rad = norm( [ state(1) state(2) state(3) ] ) ; % radius
    vel = norm( [ state(4) state(5) state(6) ] ) ; % velocity
    m = state(7) ; % mass
    
    dx = state(4) ; % velocity in x
    dy = state(5) ; % velocity in y
    dz = state(6) ; % velocity in z
    ddx = -(mu*state(1)/rad^3)+(thrust*state(4))/(1e3*m*vel) ; % acceleration in x
    ddy = -(mu*state(2)/rad^3)+(thrust*state(5))/(1e3*m*vel) ; % acceleration in y
    ddz = -(mu*state(3)/rad^3)+(thrust*state(6))/(1e3*m*vel) ; % acceleration in z
    mdot = -thrust/(g0*isp) ;

    dstate_dt = [ dx ; dy ; dz ; ddx ; ddy ; ddz ; mdot ] ;

    end