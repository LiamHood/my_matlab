%% Jason Dove & Liam Hood - HW: 6 and 7
% 5/19/19
clear all; close all; clc;
%% Givens and Constants
muearth = 398600 ; %km^3/s^2
mu = muearth ; %Satellite is Earth orbiting
%Given orbital values
rad = 6378 + 425 ;
h = sqrt( mu * rad ) ; %km^2/s
ecc = 0.0006123 ; %deg
RAAN = 86.7787 ; %deg
inc = 51.6425 ; %deg
ArgOfPerigee = 248.4036 ; %deg
TrueAnom = 0 ; %deg
Period = 2*pi/sqrt(mu)*rad^(3/2) ; %seconds
% Generate the R and V vectors
State = OrbitElem2State(h,mu,ecc,TrueAnom,RAAN,inc,ArgOfPerigee) ;
r = State( 1:3 ) ;
v = State( 4:6 ) ;
% This is the State Vector in ECI

%Omega Vector Initial
Omega = [0;0;0] ; %rad/s - This is in body frame

% Sun Vector
s = [ -1 0 0 ] ;

% Start time 
UTC = [ 2019 , 3 , 19 , 12 , 0 , 0 ] ;


%% Defining the Spacecraft Geometry and MOI
%Calculating Moments of Inertia
magdi = [ 0 0 -0.5 ] ; % magnetic dipole moment A*m^2

% Geometry means that x and y center is aligned with the center of mass.
% There is an offset for the z axis, as the sensor pulls it off the
% Z offset calc
% zoff = sensmass*1.5/( sensmass + 2*panmass + busmass ) ;
% coordinate center
% Moments of inertia of bus
Ixv = [ .679 , 2.337 , 35.326 , 300 ] ;
Iyv = [ .729 , 2.387 , 35.376 , 300 ] ;
Izv = [ .083 , .083 , .082 , .080 ] ;
Lv = [ 1 , 5 , 10 , 20 ] ;
zoffv = [ -.074 , 0 , 0 , .876 ] + .3 ;


jj = 4 ; % change jj to calculate different lengths
Ix = Ixv( jj ) ;
Iy = Iyv( jj ) ;
Iz = Izv( jj ) ;
L = Lv( jj ) ;
zoff = zoffv( 1 ) ;

I = [ Ix ; Iy ; Iz ] ;
ISpacecraft = I.*eye(3) ;
%Generating Geometry for Solar Pressure and Drag
% Only need areas along each normal axis, this is for a simplified model
bazp = .2^2  ;
bazn = bazp - 126e-6 ;
baxp = .2*.3 ;
baxn = baxp ;
bayp = baxp ;
bayn = baxn ;
layp = L*.05 ;
layn = layp ;
laxp = L*.0234 ;
laxn = laxp ;
eazp = .1^2 - 126e-6 ;
eazn = .1^2 ;
Areas = [ bazp bazn bayp bayn baxp baxn layp layn laxp layn eazp eazn ] ;

% Define Area Vector, FOR EVERY FACE (Position from COM)
COM = [0;0;zoff] ; %Center of mass is offset
%First Fine Position of center of face in the given system
%Panels (All positions are in meters
pbzp = [ 0 ; 0 ; 0 ] ;
pbzn = [ 0 ; 0 ; -.3 ] ;
pbyp = [ 0 ; .1 ; -.15 ] ;
pbyn = [ 0 ; -.1 ; -.15 ] ;
pbxp = [ .1 ; 0 ; -.15 ] ;
pbxn = [ -.1 ; 0 ; -.15 ] ;
plyp = [ 0 ; .0125 ; -.3 - L/2  ] ;
plyn = [ 0 ; -.0125 ; -.3 - L/2 ] ;
plxp = [ .0125 ; 0 ; -.3 - L/2 ] ;
plxn = [ -.0125 ; 0 ; -.3 - L/2  ] ;
pezp = [ 0 ; 0 ; -.3 - L ] ;
pezn = [ 0 ; 0 ; -.3 - L ] ;
Positions = [ pbzp pbzn pbyp pbyn pbxp pbxn plyp plyn plxp plxn pezp pezn ] - COM ;

% All faces's areas defined
%Define Normal Vectors for each
nx = [1;0;0] ;
ny = [0;1;0] ; 
nz = [0;0;1] ;
Norms = [ nz , -nz , ny , -ny , nx , -nx , ny , -ny , nx , -nx , nz , -nz ] ; %Normal Vectors for each face
% The full positions and geometries of the faces are defined

%Pressure and drag will use these areas and the rotation matrix needed to
%calculate the cross sectional area normal to the incoming force (assumed
%to be even across entire spacecraft).

%Find Initial Attitude relative to ECI
% body to LVLH
eps = [ 1 0 0 ] ;
eta = 1 ;
quat0b = [ eta , eps ] ;
Cblvlh = quat2dcm( quat0b ) ;

% LVLH to ECI
Clvlhg = ECItoLVLH( r , v ) ;
quat0 = dcm2quat( Cblvlh*Clvlhg )' ;
% Initial Euler
[Yaw0 , Pitch0, Roll0]  = quat2angle(quat0') ; %Finds initial Euler Angles
Angles = [Roll0 , Pitch0, Yaw0]' ; %Places them into vector

%% General ODE Solver - Must Contain All Perturbations
% Establish State
State = [quat0;Omega;r;v;Angles;quat0b';[0;0;0]] ;
options = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
tend = 92.75*60*10 ; 
% Set Torque
torques = [ 1 ; 1 ; 1 ; 0 ] ;
% figure(6)
% hold on
% [Time,StateNew] = ode45(@PerturbationSolver,[0,tend],State,options,I,Areas,Positions,Norms,magdi,muearth,UTC,torques,1);
% title( 'Magnitude of Torques on S/C' )
% xlabel( 'Time [Hours]' )
% ylabel( 'Torque [N-m]' )
% legend( 'Drag Torque' , 'Gravity Gradient' , 'Location' , 'eastoutside')
% % xlim([0 , 5*100/60])
% hold off
% % Set Torque
% 
figure(6)
hold on
[TimeSimple,StateNewSimple] = ode45(@PerturbationSolver,[0,tend],State,options,I,Areas,Positions,Norms,magdi,muearth,UTC,torques,0);
title( 'Magnitude of Torques on S/C' )
xlabel( 'Time [Hours]' )
ylabel( 'Torque [N-m]' )
% legend( 'Drag Torque' , 'SRP' , 'Gravity Gradient' , 'Magnetic' , 'Location' , 'eastoutside' )
legend( 'Drag Torque' , 'SRP' , 'Gravity Gradient' , 'Location' , 'eastoutside' )
% xlim([0 , 5*100/60])
hold off

%% Plots
%% Without advanced atmospheric model
% Change Time to Hours
TimeSimple = TimeSimple/3600 ; %Converts time to hours
% Relative to LVLH
% Quaternions Body to LVLH
    newquattotSimple = StateNewSimple(:,1:4) ;
    figure(1)
%     subplot(2,1,1)
    hold on
    plot(TimeSimple,StateNewSimple(:,17))
    plot(TimeSimple,StateNewSimple(:,18))
    plot(TimeSimple,StateNewSimple(:,19))
    plot(TimeSimple,StateNewSimple(:,20))
    title( 'Simple Atmosphere: Quaternions Body to LVLH' )
    legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' , 'Location' , 'eastoutside')
    xlabel('Time [Hours]')
    ylabel('Quaternion Value')
%     xlim([0 , 5*100/60]) % Limits this plot to 5 orbits
    hold off
% Euler angles Body to LVLH
    figure(2)
%     subplot(2,1,1)
    hold on
    plot(TimeSimple,rad2deg(StateNewSimple(:,21)))
    plot(TimeSimple,rad2deg(StateNewSimple(:,22)))
    plot(TimeSimple,rad2deg(StateNewSimple(:,23)))
    title( 'Simple Atmosphere: Euler Angle Body to LVLH' )
    legend( 'Roll' , 'Pitch' , 'Yaw' , 'Location' , 'eastoutside' )    
    xlabel('Time [Hours]')
    ylabel('Angle [Degrees]')
%     xlim([0 , 5*100/60]) % Limits this plot to 5 orbits
    hold off

% % Relative to ECI
% % Quaternions Body to ECI
%     figure(3)
% %     subplot(2,1,1)
%     hold on
%     plot(TimeSimple,StateNewSimple(:,1))
%     plot(TimeSimple,StateNewSimple(:,2))
%     plot(TimeSimple,StateNewSimple(:,3))
%     plot(TimeSimple,StateNewSimple(:,4))
%     title( 'Simple Atmosphere: Quaternions Body to ECI' )
%     legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' , 'Location' , 'eastoutside' )
%     xlabel('Time [Hours]')
%     ylabel('Quaternion Value')
% %     xlim([0 , 5*100/60]) % Limits this plot to 5 orbits
%     hold off
    
% % Euler Angles Body to ECI
%     figure(4)
% %     subplot(2,1,1)
%     hold on
%     plot(TimeSimple,rad2deg(StateNewSimple(:,14)))
%     plot(TimeSimple,rad2deg(StateNewSimple(:,15)))
%     plot(TimeSimple,rad2deg(StateNewSimple(:,16)))
%     title( 'Simple Atmosphere: Euler Angle Body to ECI' )
%     legend( 'Roll' , 'Pitch' , 'Yaw' , 'Location' , 'eastoutside')  
%     xlabel('Time [Hours]')
%     ylabel('Angle [Degrees]')
% %     xlim([0 , 5*100/60]) % Limits this plot to 5 orbits
%     hold off
    

% % Angular Velocity Body to ECI
%     figure(5)
% %     subplot(2,1,1)
%     for pp = 1:length(TimeSimple)
%         h_angS(:,pp) = I.*eye(3)*StateNewSimple(pp,5:7)' ; %Angular Momentums-In X, Y, and Z
%         h_MagS(pp) = norm(h_angS(:,pp)) ; % Angular Momentum Magnitude
%     end
%     plot(TimeSimple,h_angS(1,:))
%     hold on
%     plot(TimeSimple,h_angS(2,:))
%     plot(TimeSimple,h_angS(3,:))
%     plot(TimeSimple,h_MagS,'LineWidth',3)
%     title( 'Simple Atmosphere: Absolute Angular Momentum' )
%     legend( 'x-Component' , 'y-Component' , 'z-Component' , 'Magnitude' , 'Location' , 'eastoutside')
%     xlabel('Time [Hours]')
%     ylabel('Angular Momentum [kg-m^2/s]')
% %     xlim([0 , 24])
%     hold off

%% Advanced Atmospheric Model
% Time = Time/3600 ; %Converts time to hours
% % Quaternions Body to LVLH
%     figure(1)
%     subplot(2,1,2)
%     hold on
%     plot(Time,StateNew(:,17))
%     plot(Time,StateNew(:,18))
%     plot(Time,StateNew(:,19))
%     plot(Time,StateNew(:,20))
%     title( 'NRLMSISE-00 Model: Quaternions Body to LVLH' )
%     legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' , 'Location' , 'eastoutside')
%     xlabel('Time [Hours]')
%     ylabel('Quaternion Value')
%     xlim([0 , 5*100/60]) % Limits this plot to 5 orbits
%     hold off
% % Euler Angles Body to LVLH
%     figure(2)
%     subplot(2,1,2)
%     hold on
%     plot(Time,rad2deg(StateNew(:,21)))
%     plot(Time,rad2deg(StateNew(:,22)))
%     plot(Time,rad2deg(StateNew(:,23)))
%     title( 'Euler Angle LVLH' )
%     legend( 'Roll' , 'Pitch' , 'Yaw' , 'Location' , 'eastoutside' )
%     xlabel('Time [Hours]')
%     ylabel('Angle [Degrees]')
%     xlim([0 , 5*100/60]) % Limits this plot to 5 orbits
%     hold off
% 
% % Relative to ECI
% % Quaternion Body to ECI
%     newquat = StateNew(:,1:4) ;
%     figure(3)
%     subplot(2,1,2)
%     hold on
%     plot(Time,newquat(:,1))
%     plot(Time,newquat(:,2))
%     plot(Time,newquat(:,3))
%     plot(Time,newquat(:,4))
%     title( 'NRLMSISE-00 Model: Quaternions Body to ECI' )
%     legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' , 'Location' , 'eastoutside')
%     xlabel('Time [Hours]')
%     ylabel('Quaternion Value')
%     xlim([0 , 5*100/60]) % Limits this plot to 5 orbits
%     hold off
% % Euler angles Body to ECI
%     figure(4)
%     subplot(2,1,2)
%     hold on
%     plot(Time,rad2deg(StateNew(:,14)))
%     plot(Time,rad2deg(StateNew(:,15)))
%     plot(Time,rad2deg(StateNew(:,16)))
%     title( 'NRLMSISE-00 Model: Euler Angle Body to ECI' )
%     legend( 'Roll' , 'Pitch' , 'Yaw' , 'Location' , 'eastoutside' ) 
%     xlabel('Time [Hours]')
%     ylabel('Angle [Degrees]')
%     xlim([0 , 5*100/60]) % Limits this plot to 5 orbits
%     hold off
%     
% 
% % Angular Velocity
%     figure(5)
%     subplot(2,1,2)
%     for pp = 1:length(Time)
%         h_angN(:,pp) = I.*eye(3)*StateNew(pp,5:7)' ; %Angular Momentums-In X, Y, and Z
%         h_MagN(pp) = norm(h_angN(:,pp)) ; % Angular Momentum Magnitude
%     end
%     plot(Time,h_angN(1,:))
%     hold on
%     plot(Time,h_angN(2,:))
%     plot(Time,h_angN(3,:))
%     plot(Time,h_MagN,'LineWidth',3)
%     title( 'NRLMSISE-00 Model: Absolute Angular Momentum' )
%     legend( 'x-Component' , 'y-Component' , 'z-Component' , 'Magnitude' , 'Location' , 'eastoutside' )
%     xlabel('Time [Hours]')
%     ylabel('Angular Momentum [kg-m^2/s]')
%     xlim([0 , 24])
%     hold off


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
% Take in ECI state and generate LVLH vectrix. Uses this vectrix to create
% direction cosine matrix to rotate ECI to LVLH
    rd = r/norm(r) ;
    vd = v/norm(v) ;
    zLVLH = -rd ;
    yLVLH = -cross( r , v )/norm( cross( r , v ) ) ;
    xLVLH = cross( yLVLH , zLVLH ) ;
    Flvlh = [ xLVLH , yLVLH , zLVLH ] ;
    Clvlhg = Flvlh' ;
end
function ControlQuat = ComRotFind(qnow,qwant)
% Input the quaternion you want to get too and the initial quat, to get the
% control quat between them.
% ReOrder the quaternions
qconj = quatconj(qwant) ; %Finds the quaternion conjugate 
ControlQuat = quatmultiply(qconj,qnow) ; % Finds the quaternion to rotate from qnow to qwant
end

%% Perturbation Functions
% All Functions Require a State Input and output a body torque, some
% require more inputs (State(1:6) is in ECI, State(7:9) is in body frame)
% Atmospheric Drag (1)
function Torque = DragFun( Cbg , vG , rG , Areas , Positions , Norms )
% Velocity and Position are pulled in as an ECI value
vG = vG*1e3 ;
vnorm = norm(vG) ; % m/s Magnitude of velocity
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
    dotted = (Norms(:,pp)'*(vb/vnorm)) ; %Dot product of norm and sun
    if dotted > 0 
        Drag = Drag + -rho*vnorm^2*dotted*Areas(pp)*(vb/vnorm) ;
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
    rb = (Cbg*r)  ; % position vector in body frame
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
    rf = Ci2f*r ; % Ecef position
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
%% Extra Credit Atmospheric Function
function Torque = NRLMSISE( Cbg , vG , rG , Areas , Positions , Norms , UTCnow )
% Velocity and Position are pulled in as an ECI value
vnorm = norm(vG) ; % m/s Magnitude of velocity
% The normal vectors on the body frame are pulled in
% Use rotation matrix to rotate from ECI to Bodyframe
vb = Cbg*vG ; %rotates the velocity vector into the body frame
% Model of Atmosphere-NRLMSISE00

% Using the atmosnrlmsise00
lla = eci2lla(rG'*1000,UTCnow) ; % Converts eci to lat long altitude, using the utc
DaysOfYear = floor((juliandate(UTCnow)-juliandate([UTCnow(1),0,0,0,0,0]))) ;
UTseconds = (juliandate(UTCnow)-juliandate([UTCnow(1:3),0,0,0]))*24*60*60 ; %Seconds thru day
[~ , rho] = atmosnrlmsise00(lla(3) , lla(1) , lla(2) , UTCnow(1) , DaysOfYear , UTseconds , 'none') ;
rho = rho(6) ; % puls out the total mass density

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

%% General ODE Function
function dstate = PerturbationSolver(t,state,MOI,Areas,Positions,Norms,MagDi,mu,UTC,Torques,Upgrade)
% State 1 -> 4 are the quaternions, scalar first
% State 5 -> 7 is the angular velocity
% State 8 -> 10 is the position in ECI in km
% State 11 -> 13 is the veloctiy in ECI in km/s
% State 14 -> 16 is the Euler angles to ECI
% State 17 -> 20 is the quaternions to LVLH
% State 21 -> 23 is the Euler angles to LVLH

%Torques Input is a vector that allows the user to change which torques
%they want to include:
% Torques(1) =  1 includes Atmospheric Drag 
% Torques(2) =  1 includes Solar Radiation Pressure
% Torques(3) =  1 includes Gravity Gradient
% Torques(4) =  1 includes Magnetic Torque

% To see progress remove suppression
 tk = t 

s = [ -1 ; 0 ; 0 ] ;
% Find the UTC
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
% Initializing Vectors
    dstate = zeros(13,1) ;
% Set up names
    Ix = MOI(1) ;
    Iy = MOI(2) ;
    Iz = MOI(3) ;
    Ax = Areas(1) ;
    Ay = Areas(2) ;
    Az = Areas(3) ;

    quat = state( 1:4 )' ;
    w = state( 5:7 ) ;
    r = state( 8:10 ) ;
    v = state( 11:13 ) ;
    quatLVLH = state(17:20)' ;
    
% Define Rotation based off of quat
    eta = quat( 1 ) ;
    eps = quat( 2:4 )' ;
    etaLVLH = quatLVLH( 1 ) ;
    epsLVLH = quatLVLH( 2:4 )' ;
    
    Cbg = quat2dcm( quat ) ;
    Clvlhg = ECItoLVLH(r,v) ;
    wlvlh = Cbg*( cross( r , v ) / norm(r)^2 ) ;
    werr = (w - wlvlh) ; 
    
% change quaternions
    epsdot = .5 * ( eta*eye(3)*w + cross( eps , w ) ) ;
    etadot = -.5 * eps' * w ;
    dstate(1:4) = [ etadot , epsdot' ]' ;
% change quaternions relative to LVLH
    epsdote = .5 * ( etaLVLH*eye(3)*werr + cross( epsLVLH , werr ) ) ;
    etadote = -.5 * epsLVLH' * werr ;
    dstate(17:20) = [ etadote , epsdote' ]' ;     

%Simplified Orbital Solver, dState 4 -> 6 is the acceleration 
dstate(11:13) = -mu*r/norm(r)^3 ; %Position to Acceleration Function
dstate(8:10) = v ; %velocity in velocity out

    % Generate Torque
    % Baseline Torque (in Body Frame)
    Torque = [0;0;0]; %N-m - Will add to this if torques are included
    % Included Torques
    if Upgrade == 1 
        if Torques(1) == 1
                DragTorqueWu = NRLMSISE( Cbg , v , r , Areas , Positions , Norms , UTCnow) ;
                plot( tk/3600 , norm(DragTorqueWu) , 'g.' )
                Torque = Torque + DragTorqueWu ; %Torque but with the NRLMSISE-00	density	model
        end
        if Torques(2) == 1
            SRPTorqueWu = SolarPressFun( Cbg , s , Areas , Positions , Norms ) ;
            plot( tk/3600 , norm(SRPTorqueWu) , 'y.' )
            Torque = Torque + SRPTorqueWu ; %Includes the torque from SRP
        end
        if Torques(3) == 1
            GGTorqueWu = GravGradFun( Cbg , Ix , Iy , Iz , r ) ;
            plot( tk/3600 , norm(GGTorqueWu) , 'b.' )
            Torque = Torque + GGTorqueWu ; %Includes the torque from GG
        end
        if Torques(4) == 1
            MagTorqueWu = MagFun( state , Cbg , UTCnow , MagDi ) ;
            plot( tk/3600 , norm(MagTorqueWu) , 'r.' )
            Torque = Torque + MagTorqueWu ; %Includes the torque from Earth's Magnetic Field
        end
    else
        if Torques(1) == 1
                DragTorque = DragFun( Cbg , v , r, Areas , Positions , Norms ) ;
                plot( tk/3600 , norm(DragTorque) , 'g.' )
                Torque = Torque + DragTorque ; %Includes the torque from drag
        end
        if Torques(2) == 1
            SRPTorque = SolarPressFun( Cbg , s , Areas , Positions , Norms ) ;
            plot( tk/3600 , norm(SRPTorque) , 'y.' )
            Torque = Torque + SRPTorque ; %Includes the torque from SRP
        end
        if Torques(3) == 1
            GGTorque = GravGradFun( Cbg , Ix , Iy , Iz , r ) ;
            plot( tk/3600 , norm(GGTorque) , 'b.' )
            Torque = Torque + GGTorque ; %Includes the torque from GG
        end
        if Torques(4) == 1
            MagTorque = MagFun( state , Cbg , UTCnow , MagDi ) ;
            plot( tk/3600 , norm(MagTorque) , 'r.' )
            Torque = Torque + MagTorque ; %Includes the torque from Earth's Magnetic Field
        end
    end

%Use Torques to find the proper omegadot in body frame
I = MOI ;
OmegaDotx = (Torque(1)-(I(3)-I(2))*w(2)*w(3))/I(1) ; %Omegadotx Function
OmegaDoty = (Torque(2)-(I(1)-I(3))*w(1)*w(3))/I(2) ; %Omegadoty Function
OmegaDotz = (Torque(3)-(I(2)-I(1))*w(2)*w(1))/I(3) ; %Omegadotz Function
dstate(5:7) = [OmegaDotx ; OmegaDoty ; OmegaDotz ] ; %Inputs the omega dot due to given torques and motion
roll = (state(14)) ;
pitch = (state(15)) ;
% Derivative for ECI Euler
dstate(14:16) = [1,sin(roll)*tan(pitch),cos(roll)*tan(pitch);0,cos(roll),-sin(roll);0,sin(roll)*sec(pitch),cos(roll)*sec(pitch)]*w ; %Euler angle change rate ; % Euler Body to ECI

% Derivative for Body to LVLH Euler
rollLVLH = (state(21)) ;
pitchLVLH = (state(22)) ;
dstate(21:23) = [1,sin(rollLVLH)*tan(pitchLVLH),cos(rollLVLH)*tan(pitchLVLH);0,cos(rollLVLH),-sin(rollLVLH);0,sin(rollLVLH)*sec(pitchLVLH),cos(rollLVLH)*sec(pitchLVLH)]*werr ; %euler angle change rate ; % Euler Body LVLH

end




