function [ Tsrp ] = SolarPressFun( Cbg , s , areas , position , norm )
% Works at Earth. Pressure must be changed to use elsewhere. Finds the
% torque from solar radiation pressure

% Only find torque if not in eclipse
    % rotate body unit vectors to eci
    
    p = 9.08e-6 ;
    A = 0 ;
    top = zeros( 1 , 3 ) ;
    % find similarity in direction between normal face vectors and sun
    % vector
for ii = 1:length( areas )
    dir = dot( s , norm( : , ii ) ) ;
    if dir > 0 
        newarea = areas(ii)*dir ;
        A = A + newarea ;
        top = top + position(:,ii)*newarea ;
    end
end
    
    F = - p * A * s ;
    cps = sum( top ) / A ;
    cps = cps' ;
    Tsrpi = cross( cps , F ) ;
    Tsrp = Cbg*Tsrpi ;

end