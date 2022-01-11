function [ a , e , inc , RAAN , aop , ta ] = RVtoCOE( R , V )
%Calculates Classical Orbital Elements given a position and velocity vector
    % a is semi-major axis in km
    % e is eccentricity
    % inc is inclination in degrees
    % RAAN is right ascesion of node in deogrees
    % aop is argument of perigee in degrees
    % to is true anomaly in degrees
    mu = 389600; %Gravity constant of Earth in km^2 / s^2
    h = cross( R , V ); %angular momentum vector of orbit 
    k = [ 0 , 0 , 1 ];
    n = cross( k , h );
   
    %calculating semi-major axis
    a = - mu / ( 2 * (( norm(V)^2 / 2 ) - ( mu / norm(R) ))) ;
    
    %calulating eccentricity vector
    ei = ( 1 / mu ) * (((( norm(V) ^ 2 ) * ( mu / norm( R ))) * R - ( dot( R , V) ) * V )) ;
    
    %calculating eccentricity
    e = norm(ei) ;
    
    %calculating inclination
    inc = acosd( h(3) / norm(h) ); 
    
    %calculating RAAN
    RAAN = acosd( n(1) / norm(n) );
        if n(2) < 0     %Quadrant check
            RAAN = 360 - RAAN ;
        end
        
    %calculating argument of perigee    
    aop = acosd( dot( n , ei ) / ( norm(n) * norm(ei) ) );
        if ei(3) < 0    %Quadrant check
            aop = 360 - aop ;
        end
        
    %calculating true anomaly    
    ta = acosd( dot( ei , R ) / ( norm(ei) * norm(R) ));
        if dot( R , V ) < 0     %Quadrant check
            RAAN = 360 - RAAN ;
        end


end

