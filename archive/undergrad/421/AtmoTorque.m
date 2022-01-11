function [ Tatmo ] = AtmoTorque( Cbg , v , Ax , Ay , Az )
% 

    % rotate body unit vectors to eci
    nx = Cbg'*[ 1 0 0 ]' ;
    ny = Cbg'*[ 0 1 0 ]' ;
    nz = Cbg'*[ 0 0 1 ]' ;
    mnx = Cbg'*[ -1 0 0 ]' ;
    mny = Cbg'*[ 0 -1 0 ]' ;
    mnz = Cbg'*[ 0 0 -1 ]' ;
    vd = v/norm(v) ;
    
    rho = 2.5e-12 ; 
    A = 0 ;
    top = zeros( 1 , 3 ) ;
    bottom = zeros( 1 , 3 ) ;
    % find similarity in direction between normal face vectors and sun
    % vector
    dirx = dot( nx , vd ) ;
    diry = dot( ny , vd ) ;
    dirz = dot( nz , vd ) ;
    dirmx = dot( mnx , vd ) ;
    dirmy = dot( mny , vd ) ;
    dirmz = dot( mnz , vd ) ;    
    % Add to wetted area percent of area exposed
    % and add to calculations for center of pressure
    if dirx > 0 
        A = A + Ax * dirx ;
        top(1:3,1) = (nx - dirx*vd)*Ax ;
        bottom(1) = Ax ;
    end
    if dirmx > 0
        A = A + Ax * dirmx ;
        top(1:3,1) = (mnx - dirmx*vd)*Ax ;
        bottom(1) = Ax ;        
    end
    if diry > 0 
        A = A + Ay * diry ;
        top(1:3,2) = (ny - diry*vd)*Ay ;
        bottom(2) = Ay ;  
    end
    if dirmy > 0
        A = A + Ay * dirmy ;
        top(1:3,2) = (mny - dirmy*vd)*Ay ;
        bottom(2) = Ay ; 
    end
    if dirz > 0 
        A = A + Az * dirz ;
        top(1:3,3) = (nz - dirz*vd)*Az ;
        bottom(3) = Az ;
    end
    if dirmz > 0
        A = A + Az * dirmz ;
        top(1:3,3) = (mnz - dirmz*vd)*Az ;
        bottom(3) = Az ;
    end
    Cd = 2.2 ;
    F = -.5 * Cd * rho * A * (norm(v)^2) * vd ;
    cp = sum( top ) / sum( bottom ) ;
    cp = cp' ;
    Tatmoi = cross( cp , F ) ;
    Tatmo = Cbg*Tatmoi ;


end