function [ Tmag ] = MagneticTorque( Cbg , mbar , r , UTC )
% Finds magnetic torque using World Magnetic Model 2015
    dyear = decyear(UTC) ; % decimal year since 2015
    Ci2f = dcmeci2ecef( 'IAU-2000/2006' , UTC ) ; % rotation matrix for eci to ecef
    rf = Ci2f'*r ; % Ecef position
    lla = ecef2lla(rf'*1000) ; % lat long alt
    [ Bned , ~ , ~ , ~ , ~ ] = wrldmagm( lla(3) , lla(1) , lla(2) , dyear ) ; % Finds mag field in North-East-Down
    spheroid = wgs84Ellipsoid( 'm' ) ; % define spheroid of earth
    [Bf(1) , Bf(2) , Bf(3) ] = ned2ecef( Bned(1) , Bned(2) , Bned(3) , lla(1) , lla(2) , lla(3) , spheroid ) ; % mag field in ecef
    B = Ci2f'*Bf' ; % mag field in eci
    Bb = Cbg*B *1e-9 ; % mag field body frame, in Tesla
    Tmag = cross( mbar , Bb )' ; % torque in body frame
end