%% Aero 215 HW 4
% Liam Hood
% Transferring from eliptical orbit with inclination to geostationary,
% zero-inclination, orbit
clc


%% Given Values
mu = 389600; %Gravity constant of Earth in km^2 / s^2
R = [ -5096.80 , 3185.50 , 3185.50 ] ; %km
V = [ -3.95 , -6.32 , 4.74 ] ; %km/s
mRf = 42157 ; %magnitude of final position in km


%% 1 COEs of the given initial orbit

    disp( 'COEs of initial orbit' )

    %Running function to calculate COEs 
    [ a , e , inc , RAAN , aop , ta ] = RVtoCOE( R , V );

    %displaying COEs
    %Semi-major axis 
    ad = [ 'The semi-major axis is ' , num2str(a) , ' km' ] ;
    disp( ad )

    %eccentricity
    ed = [ 'The eccentricity is ' , num2str(e) ] ; 
    disp( ed )

    %Inclination
    incd = [ 'The inclination ' , num2str(inc) , ' degrees' ] ;
    disp( incd )

    %RAAN
    RAANd = [ 'The right ascenion of ascending node is ' , num2str(RAAN) , ' degrees' ] ;
    disp( RAANd )

    %Argument of perigee
    aopd = [ 'The argument of perigee is ' , num2str(aop) , ' degree' ] ;
    disp( aopd )

    %True anomaly
    tad = [ 'The true anomaly ' , num2str(ta) , ' degrees' ] ;
    disp( tad )

    %Period
    if a > 0  %Period with positive value for semi-major axis
        T = 2 * pi * (( a.^3 ) ./ ( mu )).^.5 ;
    else      %Period with negative value for semi-major axis
            T = 2 * pi * (( -a.^3 ) ./ ( mu )).^.5 ;
    end
    Td = [ 'The period is ' , num2str(T) , ' seconds' ] ;
    disp( Td )
    
    
%% 2 Simple Burn Calculations
disp( 'Simple Burns' )
    % Delta V1 circulize at apogee
    sV1f = ( 2 * ( mu / ( a * ( 1 + e ) ) - ( mu / ( 2 * ( a * ( 1 + e ) ) ))))^.5; %Velocity of circular orbit 
    sV1i = ( 2 * ( mu / ( a * ( 1 + e ) ) - ( mu / ( 2 * a )  ))) ^.5; %Velocity of original orbit at apogee
    sdV1 = sV1f - sV1i ; 
    sdV1d = [ 'The delta V of the first burn, to circulize orbit, is ' , num2str(sdV1) , ' km/s' ];
    disp( sdV1d )
    
    % Delta V2 
    sdV2 = 2 * sV1f * sind( inc / 2 ) ; %delta V to change to zero inclination
    sdV2d = [ 'The delta V of the second burn, to change to zero inclination, is ' , num2str(sdV2) , ' km/s' ];
    disp( sdV2d )
    
    % Delta V3 
    sV3f = ( 2 * ( mu / ( a * ( 1 + e ) ) - ( mu / ( mRf + ( a * ( 1 + e ) ) ))))^.5 ; %Velocity of transfer orbit at perigee
    sdV3 = sV3f - sV1f ;
    sdV3d = [ 'The delta V of the third burn, to change apogee radius to geostationary orbit, is ' , num2str(sdV3) , ' km/s' ];
    disp( sdV3d )
    
    % Delta V4
    sV4f = ( 2 * ( mu / ( mRf ) - ( mu / ( mRf * 2 ))))^.5 ; %Velocity of geostationary orbit
    sV4i = ( 2 * ( mu / ( mRf ) - ( mu / ( mRf + ( a * ( 1 + e ) ) ))))^.5 ; %Velocity of transfer orbit at apogee
    sdV4 = sV4f - sV4i;
    sdV4d = [ 'The delta V of the fourth burn, to change to geostationary orbit, is ' , num2str(sdV4) , ' km/s' ];
    disp( sdV4d )
    
    % Total Delta V
    sdV = norm(sdV1) + norm(sdV2) + norm(sdV3) + norm(sdV4);
    sdVd = [ 'The total delta V is ' , num2str(sdV) , ' km/s' ];
    disp(sdVd)
    
    
%% 3 Effecient Burn Calculations
disp( 'Effecient Burns' )
    % Delta V1 lengthening apogee
    eV1f = ( 2 * ( mu / ( a * ( 1 - e )) - ( mu / ( mRf + ( a * ( 1 - e ) ) ))))^.5 ; %Velocity of transfer orbit at perigee
    eV1i = ( 2 * ( mu / ( a * ( 1 - e ) ) - ( mu / ( 2 * a )  ))) ^.5 ; %Velocity of orignal orbit
    edV1 = eV1f - eV1i;
    edV1d = [ 'Changing apogee to geosynchronous altitude takes a delta V of ' , num2str(edV1) , ' km/s' ];
    disp(edV1d)
    
    % Delta V2 plane change and circuralization
    eV2i = ( 2 * ( mu / mRf - ( mu / ( mRf + ( a * ( 1 - e ) ) ))))^.5 ; %Velocity of transfer orbit at apogee
    eV2f = ( 2 * ( mu / mRf - mu / ( 2 * mRf ))) ^.5 ; %Velolcity of circular orbit
    edV2 = loc( eV2i , eV2f , inc ) ;
    edV2d = [ 'Plane change and circuralization takes a delta V of ' , num2str(edV2) , ' km/s' ];
    disp(edV2d)
    
    % Total delta V of effecient burn plan
    edV = norm(edV1) + norm(edV2) ;
    edVd = [ 'Changing to geostatinary orbit with effecient burns takes a delta V of ' , num2str(edV) , ' km/s' ];
    disp(edVd)
    
    
%% 4 COEs of Geostationary Orbit 
disp( 'Geostationary Orbit' )
    %State vectors for geostationary orbit
    Rf = [ mRf , 0 , 0 ] ; %km
    Rfd = [ 'Position vector for geostationary orbit in km      ' , num2str(Rf) ];
    Vf = [ 0 , eV2f , 0 ] ; %km/s
    Vfd = [ 'Velocity vector for geostationary orbit in km/s    ' , num2str(Vf) ];
    disp(Rfd)
    disp(Vfd)
    
    %COEs for Geostationary Orbit
    [ af , ef , incf , RAANf , aopf , taf ] = RVtoCOE( Rf , Vf ) ;
    
    %displaying COEs
    afd = [ 'The semi-major axis is ' , num2str(af) , ' km' ] ;
    disp( afd )

    efd = [ 'The eccentricity is ' , num2str(ef) ] ; 
    disp( efd )

    incfd = [ 'The inclination ' , num2str(incf) , ' degrees' ] ;
    disp( incfd )

    RAANfd = [ 'The right ascenion of ascending node is ' , num2str(RAANf) , ' degrees' ] ;
    disp( RAANfd )

    aopfd = [ 'The argument of perigee is ' , num2str(aopf) , ' degree' ] ;
    disp( aopfd )

    tafd = [ 'The true anomaly ' , num2str(taf) , ' degrees' ] ;
    disp( tafd )
    
    disp( 'NaN degrees means that the element is undefined due to the orbit being perfectly circular and equitorial' )

    %Calculating periods
    if a > 0
        Tf = 2 * pi * (( a.^3 ) ./ ( mu )).^.5 ;
    else
            Tf = 2 * pi * (( -a.^3 ) ./ ( mu )).^.5 ;
    end
    Tfd = [ 'The period is ' , num2str(T) , ' seconds' ] ;
    disp( Tfd )