%% Jason Dove & Liam Hood - HW: 6 and 7
% 5/19/19
clear ; close all; clc;
%% Givens and Constants
muearth = 398600 ; %km^3/s^2
mu = muearth ; %Satellite is Earth orbiting
%Given orbital values
rad = 6378 + 590 ;
h = sqrt( mu * rad ) ; %km^2/s
ecc = 0 ; %deg
RAAN = 0 ; %deg
inc = 40 ; %deg
ArgOfPerigee = 0 ; %deg
TrueAnom = 0 ; %deg
Period = 2*pi/sqrt(mu)*rad^(3/2) ; %seconds
% Generate the R and V vectors
State = OrbitElem2State(h,mu,ecc,TrueAnom,RAAN,inc,ArgOfPerigee) ;
r = State( 1:3 ) ;
v = State( 4:6 ) ;
% This is the State Vector in ECI

%Omega Vector Initial
Omega = [0;-0.001;0.002] ; %rad/s - This is in body frame

% Sun Vector
s = [ -1 0 0 ] ;

% Start time 
UTC = [ 2019 , 3 , 18 , 12 , 0 , 0 ] ;


%% Defining the Spacecraft Geometry and MOI
magdi = [ 0 0 -0.5 ] ; % magnetic dipole moment A*m^2

%Calculating Moments of Inertia
mass = 100 ;
panmass = .2*mass ; % mass of bus in kg
sensmass = .1*mass ; % sensor mass
busmass = mass - panmass - sensmass ; % solar panel mass
side = 2 ; % side length in m

% Geometry means that x and y center is aligned with the center of mass.
% There is an offset for the z axis, as the sensor pulls it off the
% Z offset calc
zoff = sensmass*1.5/(sensmass+2*panmass+busmass) ;
% coordinate center
% Moments of inertia of bus
Ixbus = (busmass*side^2)/6 + busmass*zoff^2 ;
Iybus = (busmass*side^2)/6 + busmass*zoff^2 ;
Izbus = (busmass*side^2)/6 ;

% moments of inertia of sensor
Izsens = (sensmass*.25^2)/6 ;
Iysens = (1/12)*sensmass*(.25^2+1) + sensmass*(1.5-zoff)^2 ;
Ixsens = (1/12)*sensmass*(.25^2+1) + sensmass*(1.5-zoff)^2 ;

% moments of inertia of panel
Iypan = 2*((1/12)*panmass*(.05^2+2^2) + panmass*zoff^2) ;
Izpan = 2*((1/12)*panmass*(2^2+3^2) + panmass*((2.5)^2) ) ;
Ixpan = 2*((1/12)*panmass*(.05^2+3^2) + panmass*(2.5^2+zoff^2)) ;

% moments of inertia of s/c
Ix = Ixbus + Ixsens + Ixpan ;
Iy = Iybus + Iysens + Iypan ;
Iz = Izbus + Izsens + Izpan ;
J = [ Ix ; Iy ; Iz ] ;
It = 0.6 ;
Is = 1.2 ;
mw = 1 ;
Isc = J + ( 2*It + Is )*eye(3) ;
Iw = Is*eye(3) ;

%Generating Geometry for Solar Pressure and Drag
% Only need areas along each normal axis, this is for a simplified model
APanelx = .05*3 ; %m^2 - Area of Panel x face
APanely = .05*2 ; %m^2 - Area of Panel y face
APanelz = 3*2 ; %m^2 - Area of Panel z face
ASensorx = 1*.25 ; %m^2 - Area of Sensor x face
ASensory = 1*.25 ; %m^2 - Area of Sensor y face
ASensorz = .25^2 ; %m^2 - Area of Sensor z face
ABusx = 2*2 ; %m^2 - Area of Bus faces x
ABusy = 2*2 - APanely ; %m^2 - Area of Bus face y
ABuszp = 2*2 - ASensorz ; %m^2 - Area of Bus face z positive
ABuszn = ABusx ; %m^2 - Area of Bus face z negative
Areas = [ APanelx , APanelx , APanelx , APanelx , APanely , APanely , APanelz, APanelz , APanelz , APanelz , ASensorx , ASensorx , ASensory , ASensory, ASensorz , ABusx , ABusx , ABusy , ABusy , ABuszp ,ABuszn ] ;

% Define Area Vector, FOR EVERY FACE (Position from COM)
COM = [0;0;zoff] ; %Center of mass is offset
%First Fine Position of center of face in the given system
%Panels (All positions are in meters
PPxp1 = [1;-4;0] ;
PPxp2 = [1;4;0] ;
PPxn1 = [-1;-4;0] ;
PPxn2 = [-1;4;0] ;
PPyp = [0;4;0] ;
PPyn = [0;-4;0] ;
PPzp1 = [0;-2.5;.025] ;
PPzp2 = [0;2.5;.025] ;
PPzn1 = [0;-2.5;-.025] ;
PPzn2 = [0;2.5;-.025];
PSxp = [.125;0;1.5] ;
PSxn = [-.125;0;1.5];
PSyp = [0;.125;1.5];
PSyn = [0;-.125;1.5];
PSzp = [0;0;2] ;
PBxp = [1;0;0] ;
PBxn = [-1;0;0] ;
PByp = [0;1;0] ;
PByn = [0;-1;0] ;
PBzp = [0;0;1] ;
PBzn = [0;0;-1] ;
Positions = [PPxp1 PPxp2 PPxn1 PPxn2 PPyp PPyn PPzp1 PPzp2 PPzn1 PPzn2 PSxp PSxn PSyp PSyn PSzp PBxp PBxn PByp PByn PBzp PBzn ] - COM;
% All faces's areas defined
%Define Normal Vectors for each
nx = [1;0;0] ; ny = [0;1;0] ; nz = [0;0;1] ;
%         [PPxp1  PPxp2   PPxn1   PPxn2   PPyp    PPyn    PPzp1   PPzp2   PPzn1   PPzn2   PSxp    PSxn    PSyp    PSyn    PSzp    PBxp    PBxn    PByp    PByn    PBzp    PBzn    ] - COM;

Norms = [nx     nx      -nx     -nx     ny      -ny     nz      nz      -nz     -nz     nx      -nx     ny      -ny     nz      nx      -nx     ny      -ny     nz      -nz     ] ; %Normal Vectors for each face
% The full positions and geometries of the faces are defined

%Pressure and drag will use these areas and the rotation matrix needed to
%calculate the cross sectional area normal to the incoming force (assumed
%to be even across entire spacecraft).

%Find Initial Attitude relative to ECI
% body to LVLH
eps = [ 0 0 0 ] ;
eta = 1 ;
quat0b = [ eta , eps ] ;
Cblvlh = quat2dcm( quat0b ) ;

% LVLH to ECI
Clvlhg = ECItoLVLH( r , v ) ;
quat0 = dcm2quat( Cblvlh*Clvlhg )' ;
%Initial Euler
[Yaw0 , Pitch0, Roll0]  = quat2angle(quat0') ; %Finds initial Euler Angles
Angles = [Roll0 , Pitch0, Yaw0]' ; %Places them into vector

% Control Law Constants
zeta = .7 ; % Damping coefficient
wn = .5 ;  % Wanted Natural Freq

kp = 2.*Isc.*wn^2 ;
kd = Isc.*2*zeta*wn ;
ww0 = [ 0 ; 0 ; 0 ] ;

%% General ODE Solver - Must Contain All Perturbations
% Establish State
State = [quat0;Omega;ww0;r;v;Angles;quat0b';[0;0;0];quat0] ;
options = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
tend = Period ;
% tend = 24*60*60 ;
% torques = ones( 1 , 4 ) ;
torques = [ 1 ; 1 ; 1 ; 1 ] ;
% figure
% hold on
[Time,StateNew] = ode45(@PerturbationSolver,[0,tend],State,options,Isc,Iw,Areas,Positions,Norms,magdi,muearth,UTC,torques,0,kp,kd);
% title( 'Magnitude of Torques on S/C with Upgrade Drag Model' )
% xlabel( 'Time' )
% ylabel( 'Torque (Nm)' )
% legend( 'Drag Torque' , 'SRP' , 'Gravity Gradient' , 'Magnetic' )
% hold off

%% Plots
% Relative to LVLH
    figure
    hold on
    plot(Time,StateNew(:,20))
    plot(Time,StateNew(:,21))
    plot(Time,StateNew(:,22))
    plot(Time,StateNew(:,23))
    title( 'Quaternions LVLH' )
    legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' )
    hold off
    % Euler angles
    figure
    hold on
    plot(Time,StateNew(:,24))
    plot(Time,StateNew(:,25))
    plot(Time,StateNew(:,26))
    title( 'Euler Angle LVLH' )
    legend( 'Roll' , 'Pitch' , 'Yaw' )    
    hold off

% Relative to ECI
    % LVLH to ECI
    figure
    hold on
    plot(Time,StateNew(:,27))
    plot(Time,StateNew(:,28))
    plot(Time,StateNew(:,29))
    plot(Time,StateNew(:,30))
    title( 'Quaternions LVLH to ECI' )
    legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' )
    hold off
    % Quaternion
    figure
    hold on
    plot(Time,StateNew(:,1))
    plot(Time,StateNew(:,2))
    plot(Time,StateNew(:,3))
    plot(Time,StateNew(:,4))
    title( 'Quaternions ECI' )
    legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' )
    hold off
    % Euler angles
    figure
    hold on
    plot(Time,StateNew(:,17))
    plot(Time,StateNew(:,18))
    plot(Time,StateNew(:,19))
    title( 'Euler Angle ECI' )
    legend( 'Roll' , 'Pitch' , 'Yaw' )    
    hold off
    

% Angular Velocity
figure
hold on
plot(Time,StateNew(:,5))
plot(Time,StateNew(:,6))
plot(Time,StateNew(:,7))
title( 'Angular Velocity of s/c' )
legend( 'x' , 'y' , 'z' )
hold off

% Angular Velocity
figure
hold on
plot(Time,StateNew(:,8))
plot(Time,StateNew(:,9))
plot(Time,StateNew(:,10))
title( 'Angular Velocity of wheels' )
legend( 'x' , 'y' , 'z' )
hold off

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
function Torque = SolarPressFun( Cbg , r , s , Areas , Positions , Norms )
sundirection = dot( r , s ) ;
if sundirection >= -.3 
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
else
    Torque = [ 0 ; 0 ; 0 ] ;
end
end
% Gravity Gradient (3)
function [ Tgg ] = GravGradFun( Cbg , Ix , Iy , Iz , r )
% Include all rotations to make torque into body using an inputted State
% Finds gravity gradient torque in earth orbit
    mu = 398600*(1e3)^3 ;
    rb = (Cbg*r).*1e3 ; % position vector in body frame
    rbcross = [ 0 -rb(3) rb(2) ; rb(3) 0 -rb(1) ; -rb(2) rb(1) 0 ] ;
    Ip = [ Ix , 0 , 0 ; 0 , Iy , 0 ; 0 , 0 , Iz ] ;
    Tgg = (( 3*mu )/( norm(rb)^5 ))*rbcross*Ip*rb ; % torque in body frame
end
% Earth	Magnetic Field (4)
function [ Tmag ] = MagFun( r , Cbg , UTCnow , magdi )
% Finds magnetic torque using World Magnetic Model 2015
    dyear = decyear(UTCnow) ; % decimal year since 2015
    Ci2f = dcmeci2ecef( 'IAU-2000/2006' , UTCnow ) ; % rotation matrix for eci to ecef
    rf = Ci2f'*r ; % Ecef position
    lla = ecef2lla(rf'*1000) ; % lat long alt
    [ Bned , ~ , ~ , ~ , ~ ] = wrldmagm( lla(3) , lla(1) , lla(2) , dyear ) ; % Finds mag field in North-East-Down
    spheroid = wgs84Ellipsoid( 'm' ) ; % define spheroid of earth
    [Bf(1) , Bf(2) , Bf(3) ] = ned2ecef( Bned(1) , Bned(2) , Bned(3) , lla(1) , lla(2) , lla(3) , spheroid ) ; % mag field in ecef
    B = Ci2f*Bf' ; % mag field in eci
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
vnorm = norm(vG)*1000 ; % m/s Magnitude of velocity
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
function dstate = PerturbationSolver(t,state,Isc,Iw,Areas,Positions,Norms,MagDi,mu,UTC,Torques,Upgrade,kp,kd)
% State 1 -> 4 are the quaternions, scalar first
% State 5 -> 7 is the angular velocity of s/c
% State 8 -> 10 is the angular velocity of wheels
% State 11 -> 13 is the position in ECI in km
% State 14 -> 16 is the veloctiy in ECI in km/s
% State 17 -> 19 is the Euler angles to ECI
% State 20 -> 23 is the quaternions to LVLH
% State 24 -> 26 is the Euler angles to LVLH
% State 27 -> 30 is the quaternions for LVLH to ECI

%Torques Input is a vector that allows the user to change which torques
%they want to include:
% Torques(1) =  1 includes Atmospheric Drag 
% Torques(2) =  1 includes Solar Radiation Pressure
% Torques(3) =  1 includes Gravity Gradient
% Torques(4) =  1 includes Magnetic Torque

% Remove before turning in, just for seeing progress
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

    dstate = zeros(26,1) ;
% Set up names
    Ix = Isc(1,1) ;
    Iy = Isc(2,2) ;
    Iz = Isc(3,3) ;
%     Ax = Areas(1) ;
%     Ay = Areas(2) ;
%     Az = Areas(3) ;

    quat = state( 1:4 )' ;
    wsc = state( 5:7 ) ;
    ww = state( 8:10 ) ;
    v = state( 14:16 ) ;
    r = state( 11:13 ) ;
    quatLVLH = state(20:23)' ;
    quatLVLH2eci = state(27:30)' ;
    
% Define Rotation based off of quat
    eta = quat( 1 ) ;
    eps = quat( 2:4 )' ;
    etaLVLH = quatLVLH( 1 ) ;
    epsLVLH = quatLVLH( 2:4 )' ;
    etaLVLH2eci = quatLVLH2eci( 1 ) ;
    epsLVLH2eci = quatLVLH2eci( 2:4 )' ;    
    
    Cbg = quat2dcm( quat ) ;
    Clvlhg = ECItoLVLH(r,v) ;
    cquat = quatLVLH2eci ;
    wlvlh = Clvlhg*( cross( r , v ) / norm(r)^2 ) ;
    werr = ( wsc - wlvlh ) ; 
    qerr = ComRotFind(quat,cquat) ;
    epse(1:3,1) = qerr( 2:4 ) ;
    
% change quaternions
    epsdot = .5 * ( eta*eye(3)*wsc + cross( eps , wsc ) ) ;
    etadot = -.5 * eps' * wsc ;
    dstate(1:4) = [ etadot , epsdot' ]' ;
% change quaternions relative to LVLH
    epsdote = .5 * ( etaLVLH*eye(3)*werr + cross( epsLVLH , werr ) ) ;
    etadote = -.5 * epsLVLH' * werr ;
    dstate(20:23) = [ etadote , epsdote' ]' ;   
% change quaternions LVLH to ECI
    epsdotl2e = .5 * ( etaLVLH2eci*eye(3)*wlvlh + cross( epsLVLH2eci , wlvlh ) ) ;
    etadotl2e = -.5 * epsLVLH2eci' * wlvlh ;
    dstate(27:30) = [ etadotl2e , epsdotl2e' ]' ;   

%Simplified Orbital Solver, dState 4 -> 6 is the acceleration 
dstate(14:16) = -mu*r/norm(r)^3 ; %Position to Acceleration Function
dstate(11:13) = v ; %velocity in velocity out

    % Generate Torque
    % Baseline Torque (in Body Frame)
    Torque = [0;0;0]; %N-m - Will add to this if torques are included
    % Included Torques
    if Upgrade == 1 
        if Torques(1) == 1
                DragTorqueWu = NRLMSISE( Cbg , v , r , Areas , Positions , Norms , UTCnow) ;
%                 plot( tk , norm(DragTorqueWu) , 'g.' )
                Torque = Torque + DragTorqueWu ; %Torque but with the NRLMSISE-00	density	model
        end
        if Torques(2) == 1
            SRPTorqueWu = SolarPressFun( Cbg , r , s , Areas , Positions , Norms ) ;
%             plot( tk , norm(SRPTorqueWu) , 'y.' )
            Torque = Torque + SRPTorqueWu ; %Includes the torque from SRP
        end
        if Torques(3) == 1
            GGTorqueWu = GravGradFun( Cbg , Ix , Iy , Iz , r ) ;
%             plot( tk , norm(GGTorqueWu) , 'b.' )
            Torque = Torque + GGTorqueWu ; %Includes the torque from GG
        end
        if Torques(4) == 1
            MagTorqueWu = MagFun( state , Cbg , UTCnow , MagDi ) ;
%             plot( tk , norm(MagTorqueWu) , 'r.' )
            Torque = Torque + MagTorqueWu ; %Includes the torque from Earth's Magnetic Field
        end
    else
        if Torques(1) == 1
                DragTorque = DragFun( Cbg , v , r, Areas , Positions , Norms ) ;
%                 plot( tk , norm(DragTorque) , 'g.' )
                Torque = Torque + DragTorque ; %Includes the torque from drag
        end
        if Torques(2) == 1
            SRPTorque = SolarPressFun( Cbg , r , s , Areas , Positions , Norms ) ;
%             plot( tk , norm(SRPTorque) , 'y.' )
            Torque = Torque + SRPTorque ; %Includes the torque from SRP
        end
        if Torques(3) == 1
            GGTorque = GravGradFun( Cbg , Ix , Iy , Iz , r ) ;
%             plot( tk , norm(GGTorque) , 'b.' )
            Torque = Torque + GGTorque ; %Includes the torque from GG
        end
        if Torques(4) == 1
            MagTorque = MagFun( r , Cbg , UTCnow , MagDi ) ;
%             plot( tk , norm(MagTorque) , 'r.' )
            Torque = Torque + MagTorque ; %Includes the torque from Earth's Magnetic Field
        end
    end

    dTorque = Torque ;
    cTorque = -kp*epse - kd*werr ;
    
%Use Torques to find the proper omegadot in body frame
wcross = [ 0 -wsc(3) wsc(2) ; wsc(3) 0 -wsc(1) ; -wsc(2) wsc(1) 0 ] ;
dwsc = ( Isc )\( cTorque + dTorque - wcross*Isc*wsc ) ;
dww = ( Iw )\( -cTorque - wcross*Iw*ww ) ;
dstate(5:7) = dwsc ; %Inputs the omega dot due to given torques and motion
dstate(8:10) = dww ;

roll = (state(17)) ;
pitch = (state(18)) ;
% Derivitive for ECI Euler
dstate(17:19) = [1,sin(roll)*tan(pitch),cos(roll)*tan(pitch);0,cos(roll),-sin(roll);0,sin(roll)*sec(pitch),cos(roll)*sec(pitch)]*wsc ; %Euler angle change rate ; % Euler Body to ECI

% Derivitive for Body to LVLH Euler
rollLVLH = (state(24)) ;
pitchLVLH = (state(25)) ;
dstate(24:26) = [1,sin(rollLVLH)*tan(pitchLVLH),cos(rollLVLH)*tan(pitchLVLH);0,cos(rollLVLH),-sin(rollLVLH);0,sin(rollLVLH)*sec(pitchLVLH),cos(rollLVLH)*sec(pitchLVLH)]*werr ; %uiler angle change rate ; % Euler Body LVLH

end
