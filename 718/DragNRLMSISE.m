function aDrag = DragNRLMSISE( r , v , A , m , epoch , yds)
    lla = eci2lla( r'*1e3 , epoch ) ;
    [ ~ , rho ] = atmosnrlmsise00( lla(3) , lla(1) , lla(2) , yds(1) , yds(2) , yds(3) ) ; 
    CD = 2.2 ;
    wearth = [ 0 ; 0 ; 7.29211e-5 ] ;
    vrel = v - cross( wearth , r ) ;
    aDrag = ( -.5*((CD*A)/m)*rho(6)*norm(vrel)*vrel*1e6 )*1e-3 ;
end