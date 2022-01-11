function [ tautrans , dv ] = ConstThrustIncA( ai , af , inci , incf , Isp , F , mi , steps , mu , long )
    % converts to and from canonical units internally
    alpha = [ 0 ; 2.467410607 ; -1.907470562 ; 35.892442177 ; -214.672979624 ; 947.773273608 ; -2114.861134906 ; 2271.240058672 ; -1127.457440108 ; 192.953875268 ; 8.577733773 ] ;
    beta = [ 1 ; 0.4609698838 ; 13.7756315324 ; -69.1245316678 ; 279.0671832500 ; -397.6628952136 ; -70.0139935047 ; 528.0334266841 ; -324.930383652 ; 20.583824517 ; 18.8165370778 ] ;
    
    DU = ai ; % distance unit
    TU = sqrt( ai^3 / mu ) ; % time unit sort of 
    DUoTU = sqrt( mu / ai ) ; %velocity unit
    ai = ai/DU ;
    af = af/DU ;
    mdot = -F/( Isp*9.81 ) ;
    mds = mdot/mi ;
    athrust = F/mi ;
    athrust = athrust*( 1/DU * TU^2 ) ;
    
    R = af / ai ;
    dinc = incf - inci ;
    thetaP = 0 ;
    dt = 0 ;
    dts = Period/steps ;
    
    % find theta 
    % find common node vector
    
    tol = 1e-3 ;
    aerr = 1 ;
    ierr = 1 ;
    while aerr > tol && ierr > tol
        if long == 0 
            cv = 1/( 4*lambda^2*ac^2 + 1 ) ;
        elseif long == 1 
            X = ( pi / ( 2*lambda ) )*sqrt( 1/ac ) ;
            Z = X^(-2) ;
            alphaterm = 0 ;
            betaterm = 0 ;
            for ii = 1:10 
                alphaterm = alphaterm + alpha(ii)*Z ;
                betaterm = betaterm + beta(ii)*Z ;
            end
            cv = alphaterm*betaterm ;
        end
        thetaC = atan( cos( thetaP )/ sqrt( ( 1/cv ) - 1 ) ) ;
        aNTW = athrust*[ 0 ; cos( thetaC ) ; sin( thetaC ) ] ;

        % NTW to ECI conversion
        y_NTW = v/norm( v ) ;
        z_NTW = cross( r , v )/ norm( cross( r , v ) ) ;
        x_NTW = cross( y_NTW , z_NTW ) ;
        C_ECI_NTW = [ x_NTW , y_NTW , z_NTW ] ;
        aECI = C_ECI_NTW*aNTW ;

        [ r , v ] = LowThrustCanonical( dts , r0 , v0 , mu , aECI ) ;

        % find COES

        dt = dt + dts ;
        dts = 2*pi*sqrt( acurr^3 )*( TU/steps ) ;
        
        aerr = abs( acurr - af )/af ;
        ierr = abs( inccur - incf )/incf ;
        
    end
    dvacc = 1 - sqrt( 1/R ) ;
    tf = ( -1/mds )*( 1 - exp( ( mds*dvacc )/athrusti ) ) ;
end