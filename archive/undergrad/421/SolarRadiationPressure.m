function [ Tsrp ] = SolarRadiationPressure( Cbg , s , Ax , Ay , Az )
% Works at Earth. Pressure must be changed to use elsewhere. Finds the
% torque from solar radiation pressure

% Only find torque if not in eclipse
if s ~= 0
    % rotate body unit vectors to eci
    nx = Cbg'*[ 1 0 0 ]' ;
    ny = Cbg'*[ 0 1 0 ]' ;
    nz = Cbg'*[ 0 0 1 ]' ;
    mnx = Cbg'*[ -1 0 0 ]' ;
    mny = Cbg'*[ 0 -1 0 ]' ;
    mnz = Cbg'*[ 0 0 -1 ]' ;
    
    p = 4.5e-6 ;
    A = 0 ;
    top = zeros( 1 , 3 ) ;
    bottom = zeros( 1 , 3 ) ;
    % find similarity in direction between normal face vectors and sun
    % vector
    dirx = dot( nx , s ) ;
    diry = dot( ny , s ) ;
    dirz = dot( nz , s ) ;
    dirmx = dot( mnx , s ) ;
    dirmy = dot( mny , s ) ;
    dirmz = dot( mnz , s ) ;    
    % Add to wetted area percent of area exposed
    % and add to calculations for center of pressure
    if dirx > 0 
        A = A + Ax * dirx ;
        top(1:3,1) = (nx - dirx*s)*Ax ;
        bottom(1) = Ax ;
    end
    if dirmx > 0
        A = A + Ax * dirmx ;
        top(1:3,1) = (mnx - dirmx*s)*Ax ;
        bottom(1) = Ax ;        
    end
    if diry > 0 
        A = A + Ay * diry ;
        top(1:3,2) = (ny - diry*s)*Ay ;
        bottom(2) = Ay ;  
    end
    if dirmy > 0
        A = A + Ay * dirmy ;
        top(1:3,2) = (mny - dirmy*s)*Ay ;
        bottom(2) = Ay ; 
    end
    if dirz > 0 
        A = A + Az * dirz ;
        top(1:3,3) = (nz - dirz*s)*Az ;
        bottom(3) = Az ;
    end
    if dirmz > 0
        A = A + Az * dirmz ;
        top(1:3,3) = (mnz - dirmz*s)*Az ;
        bottom(3) = Az ;
    end
    F = - p * A * s ;
    cps = sum( top ) / sum( bottom ) ;
    cps = cps' ;
    Tsrpi = cross( cps , F ) ;
    Tsrp = Cbg*Tsrpi ;
else
    Tsrp = zeros(3,1) ;
end

end