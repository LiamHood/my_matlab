close all ; clear ; clc ;
muearth = 398600 ; %km^3/s^2
mu = muearth ; %Satellite is Earth orbiting
%Given orbital values
h = 53335.2 ; %km^2/s
ecc = 0 ; %deg
RAAN = 0 ; %deg
inc = 98.43 ; %deg
ArgOfPerigee = 0 ; %deg
TrueAnom = 0 ; %deg
Period = 100*60 ; %seconds
% Generate the R and V vectors
State = OrbitElem2State(h,mu,ecc,TrueAnom,RAAN,inc,ArgOfPerigee) ;
r = State( 1:3 ) ;
v = State( 4:6 ) ;
% This is the State Vector in ECI

%Omega Vector Initial
Omega = [ 0.5 ; -2 ; 3 ]*1e-3 ; %rad/s - This is in body frame

% Sun Vector
s = [ -1 0 0 ] ;

% Start time 
UTC = [ 2019 , 3 , 19 , 12 , 0 , 0 ] ;

%% Defining the Spacecraft Geometry and MOI
%Calculating Moments of Inertia
sx = 1.5 ; % side length in m
sy = 1.5 ;
sz = 3 ;
m = 500 ; % mass of bus in kg
magdi = [ 0 0 0 ] ; % magnetic dipole moment A*m^2

% Geometry means that x and y center is aligned with the center of mass.
% There is an offset for the z axis, as the sensor pulls it off the
% Z offset calc

% moments of inertia of s/c
Ix = (1/12)*m*( sz^2 + sy^2 ) ;
Iy = (1/12)*m*( sz^2 + sx^2 ) ;
Iz = (1/12)*m*( sx^2 + sy^2 ) ;
I = [ Ix 0 0 ; 0 Iy 0 ; 0 0 Iz ] ;
%Generating Geometry for Solar Pressure and Drag
% Only need areas along each normal axis, this is for a simplified model
Areax = 3*1.5 ;
Areay = 3*1.5 ;
Areaz = 1.5*1.5 ;
Areas = [ Areax , Areax , Areay , Areay , Areaz , Areaz ] ;

% Define Area Vector, FOR EVERY FACE (Position from COM)
%First Fine Position of center of face in the given system
%Panels (All positions are in meters
Px = [ .75 , 0 , 0 ]' ;
Py = [ 0 , .75 , 0 ]' ;
Pz = [ 0 , 0 , 1.5 ]' ;
Positions = [ Px -Px Py -Py Pz -Pz ] ;
%Define Normal Vectors for each
nx = [1;0;0] ; ny = [0;1;0] ; nz = [0;0;1] ;
Norms = [ nx , -nx , ny , -ny , nz , -nz ] ; %Normal Vectors for each face
% The full positions and geometries of the faces are defined

%Pressure and drag will use these areas and the rotation matrix needed to
%calculate the cross sectional area normal to the incoming force (assumed
%to be even across entire spacecraft).

%Find Initial Attitude relative to ECI
% body to LVLH
eps = [ 0 0 .0 ] ;
eta = 1 ;
quat0b = [ eta , eps ] ;
Cblvlh = quat2dcm( quat0b ) ;

% LVLH to ECI
Clvlhg = ECItoLVLH( r , v ) ;
quat0 = dcm2quat( Cblvlh'*Clvlhg )' ;

%Initial Euler
[Yaw0 , Pitch0, Roll0]  = quat2angle(quat0') ; %Finds initial Euler Angles
Angles = [Yaw0 , Pitch0, Roll0]' ; %Places them into vector

% Control Law Constants
zeta = .7 ; % Damping coefficient
wn = .5 ;

kp = 2.*I.*wn^2 ;
kd = I.*2*zeta*wn ;
%% General ODE Solver - Must Contain All Perturbations
% Establish State
State = [quat0;Omega;r;v;Angles] ;
options = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
tend = 100*60 ;
[Time,StateNew] = ode45(@ControlAndDisturbances,[0,tend],State,options,I,Areas,Positions,Norms,magdi,muearth,UTC,[0,0,0,0],0, kp , kd );

%% Plots
for ii = 1:length( Time )
    r = StateNew( ii,8:10 )' ;
    v = StateNew( ii,11:13 )' ;
    quat = StateNew( ii,1:4 ) ;
    Clvlhg(:,:,ii) = ECItoLVLH( r , v ) ;
    cquat(ii,:) = dcm2quat( Clvlhg(:,:,ii) ) ;
    qstar = quatconj( cquat(ii,:) ) ;
    qerr(ii,:) = quatmultiply( qstar , quat ) ;
end

% Relative to ECI
    % Quaternion
    newquattot = StateNew(:,1:4) ;
    figure
    hold on
    plot(Time,newquattot(:,1))
    plot(Time,newquattot(:,2))
    plot(Time,newquattot(:,3))
    plot(Time,newquattot(:,4))
    title( 'Quaternions ECI' )
    legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' )
    hold off
    
    % Euler angles
    figure
    hold on
    plot(Time,StateNew(:,13))
    plot(Time,StateNew(:,14))
    plot(Time,StateNew(:,15))
    title( 'Euler Angle ECI' )
    legend( 'Roll' , 'Pitch' , 'Yaw' )    
    hold off
    
% Angular Velocity
figure
hold on
plot(Time,StateNew(:,5))
plot(Time,StateNew(:,6))
plot(Time,StateNew(:,7))
title( 'Angular Velocity' )
legend( 'x' , 'y' , 'z' )
hold off

    figure
    hold on
    plot(Time,qerr(:,1))
    plot(Time,qerr(:,2))
    plot(Time,qerr(:,3))
    plot(Time,qerr(:,4))
    title( 'Quaternion Error' )
    legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' )
    hold off
    
%     figure
%     hold on
%     plot(Time,cquat(:,1))
%     plot(Time,cquat(:,2))
%     plot(Time,cquat(:,3))
%     plot(Time,cquat(:,4))
%     title( 'Command Quaternions' )
%     legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' )
%     hold off
    
%     figure
%     hold on
%     for ii = 1:3
%         for jj = 1:3 
%             rot(:) = Clvlhg(ii,jj,:) ;
%             plot( Time , rot )
%         end
%     end
%     legend( '1,1' , '1,2' , '1,3' , '2,1' , '2,2' , '2,3' , '3,1' , '3,2' , '3,3' )
%     hold off
    
%     r = StateNew( :,8:10 )' ;
%     v = StateNew( :,11:13 )' ;
%     figure 
%     plot3( r(1,:) , r(2,:) , r(3,:) )
%     figure 
%     plot3( v(1,:) , v(2,:) , v(3,:) )
    
%% Helpful Functions
%Rotation Fcn in degrees
function Mat = xRot(Angle)
Mat = [1,0,0;0,cosd(Angle),sind(Angle);0,-sind(Angle),cosd(Angle)] ; %defines the x rotation matrix
end
function Mat = yRot(Angle)
Mat = [cosd(Angle),0,-sind(Angle);0,1,0;sind(Angle),0,cosd(Angle)] ; %defines the y rotation matrix
end
function Mat = zRot(Angle)
Mat = [cosd(Angle),sind(Angle),0;-sind(Angle),cosd(Angle),0;0,0,1] ; %defines the z rotation matrix
end
function State = OrbitElem2State(h,mu,ecc,TA,RAAN,inc,ArgOfPerigee)
%Input all orbital elements to get position and velocity states 
%State(1:3) = position (km)
%State(4:6) = velocity (km/s)
State(1:3,1) = h^2/(mu*(1+ecc*cosd(TA)))*[cosd(TA);sind(TA);0] ; %position in perifocal
State(4:6,1) = mu/h*[-sind(TA);(ecc+cosd(TA));0] ; %Velocity in Perifocal
%state is now in perifocal, must rotate accordingly
State(1:3,1) = zRot(-RAAN)*xRot(-inc)*zRot(-ArgOfPerigee)*State(1:3) ; %transforms the position into perifocal
State(4:6,1) = zRot(-RAAN)*xRot(-inc)*zRot(-ArgOfPerigee)*State(4:6) ; %transforms the velocity into perifocal
end

%% Reference Frame Functions
function Clvlhg = ECItoLVLH( r , v )
% Take in ECI state and generate LVLH vectrix. Use this vectrix to create
% direction cosine matrix to rotate ECI to LVLH
    rd = r/norm(r) ;
    vd = v/norm(v) ;
    zLVLH = -rd ;
    yLVLH = -cross( rd , vd ) ;
    xLVLH = cross( yLVLH , zLVLH ) ;
    Flvlh = [ xLVLH , yLVLH , zLVLH ] ;
    Clvlhg = Flvlh ;
end

%% Perturbation Functions
% All Functions Require a State Input and output a body torque, some
% require more inputs (State(1:6) is in ECI, State(7:9) is in body frame)
% Atmospheric Drag (1)
function Torque = DragFun( Cbg , vG , rG , Areas , Positions , Norms )
% Velocity and Position are pulled in as an ECI value
vnorm = norm(vG)*1000 ; % m/s Magnitude of velocity
alt = (norm(rG) - 6378) ; %km - altitude of spacecraft
% The normal vectors on the body frame are pulled in
% Use rotation matrix to rotate from ECI to Bodyframe
vb = Cbg*vG ; %rotates the velocity vector into the body frame
% Model of Atmosphere
rho = EarthStdAtm(alt) ; %inputs the altitude of the sc to get density of air
% Find the amount of area on each side, exposed to the sun, use dot products
Drag = 0 ; %preallocates force summation
Top = [0;0;0] ; %Preallocates Cps Parts 
Bot = 0 ; 
    for pp = 1:length(Areas)
        dotted = (Norms(:,pp)'*vb) ; %Dot product of norm and sun
        if dotted > 0 
            Drag = Drag + -rho*vnorm^2*vb*dotted*Areas(pp) ;
        end 
        % Need Center of Pressure Now
        if dotted > 0 
            Top = Top + Positions(:,pp) * dotted * Areas(pp) ; % Sums the position vector times the dotted value and area
            Bot = Bot + dotted * Areas(pp) ; %Sums the n dot s * area
        end
    end
Cps = Top/Bot ; %Defines the center of pressure
Torque = cross(Cps,Drag) ; %Cross product of force and distance

end
% Solar	Pressure (2)
function Torque = SolarPressFun( Cbg , s , Areas , Positions , Norms )
%At 1 AU, solar pressure is 
ps = 9.08E-6 ; %Pascals
%The normal vectors on the body frame are pulled in
%Use rotation matrix to rotate from ECI to Bodyframe
sb = Cbg*s ; %rotates the sun vector into the body frame
%Find the amount of area on each side, exposed to the sun, use dot products
Forces = 0 ; %preallocates force summation
Top = [0;0;0] ; %Preallocates Cps Parts 
Bot = 0 ; 
for pp = 1:length(Areas)
    dotted = (Norms(:,pp)'*s) ; %Dot product of norm and sun
    if dotted > 0 
        Forces = Forces + -ps*sb*dotted*Areas(pp) ;
    end 
    % Need Center of Pressure Now
    
    if dotted > 0 
        Top = Top + Positions(:,pp) * dotted * Areas(pp) ; % Sums the position vector times the dotted value and area
        Bot = Bot + dotted * Areas(pp) ; %Sums the n dot s * area
    end
end
Cps = Top/Bot ; %Defines the center of pressure
Torque = cross(Cps,Forces) ; %Cross product of force and distance
end
% Gravity Gradient (3)
function [ Tgg ] = GravGradFun( Cbg , Ix , Iy , Iz , r )
% Include all rotations to make torque into body using an inputted State
% Finds gravity gradient torque in earth orbit
    mu = 398600 ;
    rb = Cbg*r ; % position vector in body frame
    rbcross = [ 0 -rb(3) rb(2) ; rb(3) 0 -rb(1) ; -rb(2) rb(1) 0 ] ;
    Ip = [ Ix , 0 , 0 ; 0 , Iy , 0 ; 0 , 0 , Iz ] ;
    Tgg = (( 3*mu )/( norm(rb)^5 ))*rbcross*Ip*rb ; % torque in body frame
end
% Earth	Magnetic Field (4)
function [ Tmag ] = MagFun( state , Cbg , UTCnow , magdi )
% Finds magnetic torque using World Magnetic Model 2015
r = state( 8:10 ) ;
    dyear = decyear(UTCnow) ; % decimal year since 2015
    Ci2f = dcmeci2ecef( 'IAU-2000/2006' , UTCnow ) ; % rotation matrix for eci to ecef
    rf = Ci2f'*r ; % Ecef position
    lla = ecef2lla(rf'*1000) ; % lat long alt
    [ Bned , ~ , ~ , ~ , ~ ] = wrldmagm( lla(3) , lla(1) , lla(2) , dyear ) ; % Finds mag field in North-East-Down
    spheroid = wgs84Ellipsoid( 'm' ) ; % define spheroid of earth
    [Bf(1) , Bf(2) , Bf(3) ] = ned2ecef( Bned(1) , Bned(2) , Bned(3) , lla(1) , lla(2) , lla(3) , spheroid ) ; % mag field in ecef
    B = Ci2f'*Bf' ; % mag field in eci
    Bb = Cbg*B *1e-9 ; % mag field body frame, in Tesla
    Tmag = cross( magdi , Bb )' ; % torque in body frame
end


%% Earth Standard Atmosphere - Chart with linear interpolation
function [rho] = EarthStdAtm(alt)
% Chart Values Provided by the website below
% http://toughsf.blogspot.com/2017/09/low-earth-orbit-atmospheric-scoops.html
Altitudes = [50 60 70 80 100 200 300 400 500 600 700 800 900 1000] ; % Chart altitudes in km
Densities = [1.03E-3 3.10E-4 8.28E-5 1.85E-5 5.36E-7 3.13E-10 2.40E-11 3.38E-12 6.21E-13 1.39E-13 4.03E-14 1.66E-14 9.11E-15 5.85E-15] ; %Densities for given altitudes in kg/m^3
% Use linear interpolation
Pos = find(Altitudes>alt) ;
if alt>1000
    warning('AHHHHHHH, altitude too high!')
end
Dens1 = Densities(Pos(1)) ; %Finds the density that is right below input
Dens2 = Densities(Pos(1)+1) ; %Finds density above input altitude
AltDiff = Altitudes(Pos(1)+1)-Altitudes(Pos(1)) ; % Finds Difference in altitudes
rho = Dens1 + (Dens2-Dens1)/AltDiff ; %Linearly interpolates the rho at a given altitude

end

%% General ODE Function
function dstate = ControlAndDisturbances(t,state,MOI,Areas,Positions,Norms,MagDi,mu,UTC,Torques,Upgrade,kp,kd)
% State 1 -> 4 are the quaternions, scalar first
% State 5 -> 7 is the angular velocity
% State 8 -> 10 is the position in ECI in km
% State 11 -> 13 is the veloctiy in ECI in km/s
%Torques Input is a vector that allows the user to change which torques
%they want to include:
% Torques(1) =  1 includes Atmospheric Drag 
% Torques(2) =  1 includes Solar Radiation Pressure
% Torques(3) =  1 includes Gravity Gradient
% Torques(4) =  1 includes Magnetic Torque

% Remove before turning in, just for seeing progress
t

s = [ -1 ; 0 ; 0 ] ;
    
     if t >= 60
         tmin = floor( t/60 ) ;
         t = t - tmin*60 ;
         if tmin >= 60
             thour = floor( tmin/60 ) ;
             tmin = tmin - thour*60 ;
         else
             thour = 0 ;
         end
     else
         tmin = 0 ;
         thour = 0 ;
     end
     if ( thour + UTC(4) ) >= 24 
         thour = thour - 24 ;
         UTC(3) = UTC(3) + 1 ;
     end
     UTCnow = UTC + [ 0 0 0 thour tmin t ] ;
    JDay = juliandate( UTCnow ) ; % Find Julian day of current moment
    J2000 = JDay - juliandate( [ 2000 , 1 , 1 , 12 , 0 , 0 ] ) ;

    dstate = zeros(13,1) ;
% Set up names
    Ix = MOI(1,1) ;
    Iy = MOI(2,2) ;
    Iz = MOI(3,3) ;


    quat = state( 1:4 )' ;
    w = state( 5:7 ) ;
    r = state( 8:10 ) ;
    v = state( 11:13 ) ;
    
    eta = quat( 1 ) ;
    eps(:,1) = quat( 2:4 ) ;
    Cbg = quat2dcm( quat ) ;

    Clvlhg = ECItoLVLH( r , v ) ;
    cquat = dcm2quat( Clvlhg' ) ;
    quat = [ eta ; eps ]' ;
    qstar = quatconj( cquat ) ;
    qerr = quatmultiply( qstar , quat ) ;
    epse(1:3,1) = qerr( 2:4 ) ;
    wlvlh = Clvlhg*( cross( r , v ) / norm(r)^2 ) ;
    werror = w - wlvlh ; 
    
% change quaternions
    epsdot = .5 * ( eta*eye(3)*w + cross( eps , w ) ) ;
    etadot = -.5 * eps' * w ;
    dstate(1:4) = [ etadot , epsdot' ]' ;
        

%Simplified Orbital Solver, dState 4 -> 6 is the acceleration 
dstate(11:13) = -mu*r/norm(r)^3 ; %Position to Acceleration Function
dstate(8:10) = v ; %velocity in velocity out

    % Generate Torque
    % Baseline Torque (in Body Frame)
    Torque = [0;0;0]; %N-m - Will add to this if torques are included
    % Included Torques
    if Torques(1) == 1
        if Upgrade == 1 %if upgrade is 1, then the new atmospheric model is improved
            Torque = Torque + NRLMSISEFun(State,Cbg,Areas) ; %Torque but with the NRLMSISE-00	density	model
        else
            Torque = Torque + DragFun( Cbg , v , r, Areas , Positions , Norms ) ; %Includes the torque from drag
        end
    end
    if Torques(2) == 1
        Torque = Torque + SolarPressFun( Cbg , s , Areas , Positions , Norms ) ; %Includes the torque from SRP
    end
    if Torques(3) == 1
        Torque = Torque + GravGradFun( Cbg , Ix , Iy , Iz , r ) ; %Includes the torque from GG
    end
    if Torques(4) == 1
        if norm( MagDi ) ~= 0
            Torque = Torque + MagFun( state , Cbg , UTCnow , MagDi ) ; %Includes the torque from Earth's Magnetic Field
        else
            Torque = Torque ;
        end
    end
    CTorque = -kp*epse - kd*werror ; 
    
%Use Torques to find the proper omegadot in body frame
I = MOI ;
wcross = [ 0 -w(3) w(2) ; w(3) 0 -w(1) ; -w(2) w(1) 0 ] ;
dw = inv( I )*( -wcross*I*w + Torque + CTorque ) ;
dstate(5:7) = dw ; %Inputs the omega dot due to given torques and motion

dstate(14:16) = state(5:7) ;

    function Clvlhg = ECItoLVLH( r , v )
    % Take in ECI state and generate LVLH vectrix. Use this vectrix to create
    % direction cosine matrix to rotate ECI to LVLH
        rd = r/norm(r) ;
        vd = v/norm(v) ;
        zLVLH = -rd ;
        yLVLH = -cross( r , v )/norm( cross( r , v ) ) ;
        xLVLH = cross( yLVLH , zLVLH ) ;
        Flvlh = [ xLVLH , yLVLH , zLVLH ] ;
        Clvlhg = Flvlh' ;
    end

end