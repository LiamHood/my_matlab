A%% HW 9 -- Jason Dove - Liam Hood
%6/7/19
%% Clear Up
close all ; clear all ; clc ;

%% Givens and Knowns
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
StateOrb = OrbitElem2State(h,mu,ecc,TrueAnom,RAAN,inc,ArgOfPerigee) ;
r = StateOrb( 1:3 ) ;
v = StateOrb( 4:6 ) ;
% This is the State Vector in ECI

%Omega Vector Initial
wsc0 = [ 0.5 ; -2 ; 3 ]*1e-3 ; %rad/s - This is in body frame
ww0 = [ 0 ; 0 ; 0 ] ;

%% Defining the Spacecraft Geometry and MOI
%Calculating Moments of Inertia
sx = 1.5 ; % side length in m
sy = 1.5 ; 
sz = 3 ;
m = 500 ; % mass of bus in kg

% Rxn Wheel info
It = 0.6 ;
Is = 1.2 ;
mw = 1 ;

% moments of inertia of s/c
Ix = (1/12)*m*( sz^2 + sy^2 ) ;
Iy = (1/12)*m*( sz^2 + sx^2 ) ;
Iz = (1/12)*m*( sx^2 + sy^2 ) ;
J = [ Ix 0 0 ; 0 Iy 0 ; 0 0 Iz ] ;
Isc = J + ( 2*It + Is )*eye(3) ;
Iw = Is*eye(3) ;

%Find Initial Attitude relative to ECI
% body to LVLH
eps = [ 0 0 0 ] ;
eta = 1 ;
quat0b = [ eta , eps ] ;
Cblvlh = quat2rotm( quat0b ) ;

% LVLH to ECI
Clvlhg = ECItoLVLH( r , v ) ;
% Quat defining ECI to Body
quat0 = dcm2quat( Cblvlh*Clvlhg )' ;

%Initial Euler
[Yaw0 , Pitch0, Roll0]  = quat2angle( quat0' ) ; %Finds initial Euler Angles
Angles = [Roll0 , Pitch0, Yaw0]' ; %Places them into vector

% Control Law Constants
zeta = .7 ; % Damping coefficient
wn = .5 ;  % Wanted Natural Freq

% Found Gains
kp = 2.*Isc.*wn^2 ;
kd = Isc.*2*zeta*wn ;

%% General ODE Solver - Must Contain All Perturbations
% Establish State
State = [ quat0 ; wsc0 ; ww0 ; r ; v ; Angles ; quat0 ] ;
options = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
tend = 60*100 ;
[ Time , StateNew ] = ode45( @RxnWheels , [0,tend] , State , options , Isc , Iw , muearth , kp , kd ) ;
% State 1 -> 4 are the quaternions, scalar first
% State 5 -> 7 is the angular velocity of s/c
% State 8 -> 10 is the angular velocity of wheels
% State 11 -> 13 is the position in ECI in km
% State 14 -> 16 is the veloctiy in ECI in km/s
% State 17 -> 19 is the Euler Angles in ECI in 


%% Plots
%% Plots
Time = Time/(60^2) ;
for ii = 1:length( Time )
    hsc(:,ii) = Isc*StateNew(ii,5:7)' ;
end
% Relative to LVLH
% Quaternion
    figure
    subplot( 2 , 1 , 1 )
    hold on
    plot(Time,StateNew(:,20))
    plot(Time,StateNew(:,21))
    plot(Time,StateNew(:,22))
    plot(Time,StateNew(:,23))
    title( 'Quaternions of LVLH Relative to ECI' )
    legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' , 'Location' , 'eastoutside')
    ylabel( 'Quaternion Value Magnitude' )
    xlabel( 'Time [Hours]' )
    hold off
    
% Relative to ECI
    % Quaternion
    subplot( 2 , 1 , 2 )
    hold on
    plot(Time,StateNew(:,1))
    plot(Time,StateNew(:,2))
    plot(Time,StateNew(:,3))
    plot(Time,StateNew(:,4))
    title( 'Quaternions of Spacecraft Relative to ECI' )
    legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' , 'Location' , 'eastoutside' )
    ylabel( 'Quaternion Value Magnitude' )
    xlabel( 'Time [Hours]' )
    hold off
    
    % Euler angles
    figure
    subplot(2,1,1)
    hold on
    plot(Time,StateNew(:,17)*(180/pi))
    plot(Time,StateNew(:,18)*(180/pi))
    plot(Time,StateNew(:,19)*(180/pi))
    title( 'Euler Angle of Spacecraft Relative to ECI' )
    legend( 'Roll' , 'Pitch' , 'Yaw' ,'Location' , 'eastoutside')  
    ylabel( 'Euler Angles [Degrees]' )
    xlabel( 'Time [Hours]' )    
    hold off
    
    
    % Angular Velocity of S/C
    StateNew(ii,5:7)
    subplot(2,1,2)
    hold on
    plot(Time,StateNew(:,5)*(180/pi))
    plot(Time,StateNew(:,6)*(180/pi))
    plot(Time,StateNew(:,7)*(180/pi))
    title( 'Angular Velocity of Spacecraft' )
    legend( 'x-axis' , 'y-axis' , 'z-axis' ,'Location' , 'eastoutside')
    ylabel( 'Angular Momentum [deg/s]' )
    xlabel( 'Time [Hours]' )   
    hold off

    % Angular Momentum of S/C
    for pp = 1:length(Time)
        normalhsc(pp) = norm(hsc(:,pp)) ;
    end
    figure
    hold on
    plot(Time,hsc(1,:))
    plot(Time,hsc(2,:))
    plot(Time,hsc(3,:))
    plot(Time,normalhsc,'LineWidth',3)
    title( 'Angular Momentum of Spacecraft' )
    legend( 'x-axis' , 'y-axis' , 'z-axis' , 'Magnitude' )
    ylabel( 'Angular Momentum [kg*m^2/s]' )
    xlabel( 'Time [Hours]' )   
    hold off

    % Angular Velocity of Reaction Wheels
    figure
    hold on
    plot(Time,StateNew(:,8)*(1/(2*pi)))
    plot(Time,StateNew(:,9)*(1/(2*pi)))
    plot(Time,StateNew(:,10)*(1/(2*pi)))
    title( 'Spin Rate of Each Axis'' Reaction Wheel' )
    legend( 'x-axis' , 'y-axis' , 'z-axis' )
    ylabel( 'Spin Rate [revolutions/s]' )
    xlabel( 'Time [Hours]' )
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
%% Command Rotation Finder
function ControlQuat = ComRotFind(qnow,qwant)
% Input the quaternion you want to get too and the initial quat, to get the
% control quat between them.
% ReOrder the quaternions
qconj = quatconj(qwant) ; %Finds the quaternion conjugate 
ControlQuat = quatmultiply(qconj,qnow) ; % Finds the quaternion to rotate from qnow to qwant
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
%% HomeMade Quat2Rotm
function C = quat2rotm(quat)
% Input the quat in the [eta , eps ] form
eta = quat(1) ; 
eps = quat(2:4)' ;
epsx = [0,-eps(3),eps(2);eps(3),0,-eps(1);-eps(2),eps(1),0] ;
C = (2*eta-1)*eye(3)+2*eps*eps'-2*eta*epsx ; % Defines the rotation matrix for the given quat
end
%% HomeMade Rotm2Quat
function quat = rotm2quat(C)
% Input 3x3 rotation and get [eta , eps] out
eta = sqrt(trace(C)+1)/2 ;
eps(1) = (C(2,3)-C(3,2))/(4*eta) ;
eps(2) = (C(3,1)-C(1,3))/(4*eta) ;
eps(3) = (C(1,2)-C(2,1))/(4*eta) ;
quat = [eta , eps] ;
end
%% General ODE Function
function dstate = RxnWheels(t,state,Isc,Iw,mu,kp,kd)
% State 1 -> 4 are the quaternions, scalar first
% State 5 -> 7 is the angular velocity of s/c
% State 8 -> 10 is the angular velocity of wheels
% State 11 -> 13 is the position in ECI in km
% State 14 -> 16 is the veloctiy in ECI in km/s
% State 17 -> 19 is the Euler Angles in ECI in 
% State 20 -> 23 is the quaternion of ECI to LVLH
% Remove before turning in, just for seeing progress
t

    dstate = zeros(13,1) ;
% Set up names

    quat = state( 1:4 )' ;
    wsc = state( 5:7 ) ;
    ww = state( 8:10 ) ;
    r = state( 11:13 ) ;
    v = state( 14:16 ) ;
    quatLVLH = state(20:23)' ;
    
% Define Rotation based off of quat
    eta = quat( 1 ) ;
    eps = quat( 2:4 )' ;
    etaLVLH = quatLVLH( 1 ) ;
    epsLVLH = quatLVLH( 2:4 )' ;
% Determine LVLH frame
    Clvlhg = ECItoLVLH( r , v ) ;  
    cquat = quatLVLH  ;
    CbG = quat2dcm(quat) ;
    qerr = ComRotFind(quat,cquat) ;
    epse(1:3,1) = qerr( 2:4 ) ;
    wlvlh = CbG*( cross( r , v ) / norm(r)^2 ) ;

    werror = (wsc - wlvlh) ; 
    
% change quaternions body to eci
    epsdot = .5*(eta*eye(3)+[0,-eps(3),eps(2);eps(3),0,-eps(1);-eps(2),eps(1),0])*wsc ;
    etadot = -.5*eps'*wsc ;
    dstate(1:4) = [ etadot , epsdot' ]' ;
% eci to LVLH
    epsdotLVLH = .5*(etaLVLH*eye(3)+[0,-epsLVLH(3),epsLVLH(2);epsLVLH(3),0,-epsLVLH(1);-epsLVLH(2),epsLVLH(1),0])*wlvlh ;
    etadotLVLH = -.5*epsLVLH'*wlvlh ;
    dstate(20:23) = [etadotLVLH , epsdotLVLH']  ;

%Simplified Orbital Solver, dState 4 -> 6 is the acceleration 
dstate(14:16) = -mu*r/norm(r)^3 ; %Position to Acceleration Function
dstate(11:13) = v ; %velocity in velocity out

% Define the control Torque needed
    %CTorque = 0;
    CTorque = -kp*epse - kd*werror ;
%Use Torques to find the proper omegadot in body frame
wcross = [ 0 -wsc(3) wsc(2) ; wsc(3) 0 -wsc(1) ; -wsc(2) wsc(1) 0 ] ;
dwsc = ( Isc )\( CTorque - wcross*Isc*wsc ) ;
dww = ( Iw )\( -CTorque - wcross*Iw*ww ) ;
dstate(5:7) = dwsc ; %Inputs the omega dot due to given torques and motion
dstate(8:10) = dww ;

% Euler Angles Change Rates
roll = state(17) ;
pitch = state(18) ;
dstate(17:19) = [1,sin(roll)*tan(pitch),cos(roll)*tan(pitch);0,cos(roll),-sin(roll);0,sin(roll)*sec(pitch),cos(roll)*sec(pitch)]*wsc ; %Euler angle change rate ; % Euler Body to ECI

end


