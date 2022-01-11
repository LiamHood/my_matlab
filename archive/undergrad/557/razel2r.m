function [ r ] = razel2r( rho , az , el , lat , long , alt , utc )
    rsite_ecef = lla2ecef( [ lat , long , alt ] )'*1e-3 ;
    rhosez = [ -rho*cos(el)*cos(az) ; rho*cos(el)*sin(az) ; rho*sin(el) ] ;
    sez2ecef = [ sind(lat)*cosd(long) , -sind(long) , cosd(lat)*cosd(long) ; ...
                sind(lat)*cosd(long) , cosd(long) , cosd(lat)*sind(long) ; ...
                -cosd(lat) , 0 , sind(lat) ] ;
    rhoecef = sez2ecef*rhosez ;
    eci2ecef = dcmeci2ecef( 'IAU-2000/2006' , utc ) ;
    recef = rsite_ecef + rhoecef ;
    r = eci2ecef'*recef ;
end