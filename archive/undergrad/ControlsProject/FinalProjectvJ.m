%% Jason Dove & Liam Hood - Final Project
% 5/19/19
clear ; close all; clc;
%% Givens and Constants
muearth = 398600 ; %km^3/s^2
mu = muearth ; %Satellite is Earth orbiting
%Given orbital values
rad = 6378 + 400 ;
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

%Omega Vector of Detumble    
OmegaDet = [0.05 ; -0.02 ; -0.03]*1e1 ; %rad/s - Body Frame

%Omega Vector Initial (For Beginning of Operational)
Omega = [0;-0.001;0.002] ; %rad/s - This is in body frame

% Sun Vector
s = [ -1 0 0 ] ;

% Start time 
UTC = [ 2019 , 3 , 18 , 12 , 0 , 0 ] ;


%% Attitude Determination
% Assignment givens were lacking the sensor readouts, so we assumed ours
% Use Quest Method
% GPS Accuracy: https://www.hindawi.com/journals/mpe/2012/596396/
%       "main error sources affecting attitude determination are the
        % multipath, receiver internal noise, and antenna phase center variation "
        % GPS is much less accurate than the Star Trackers, so we set the weighting
        % factor to w = 1/sigma^2 (sigma is assumed to be accuracy)
        
        %Assume can get 1 degree of accuracy or 0.017453292519943 radians
% 
% Star Tracker Accuracy:
%       https://engineering.purdue.edu/~bethel/ct601.pdf
%       "Angular accuracy
        % (Full temperature range), per star
        % coordinate, per update (arc sec, 1 ?)  = 3
        % Noise equivalant angleb (arc sec, 1 ?) = 5
        
        %Total Error of 4.84814e-6*8 radians =  3.878512000000000e-05
        %radians
        
        % Star Tracker is 400 times more accuracte than GPS 
        % We will set the weight factor to w = 1/sigma^2 , where sigma is
        % the error (1 sigma)
        
% We will place 2 star trackers, one on each side without solar panels, but
% not across from the. Assume the Star Trackers will give us the position
% of the point of Aries and the the Z-axis of the ECI frame to the spacecraft. We know that
% a star tracker can fully determine attitude, because it can take the
% orientation of a starfield and back out total attitude from that.

% We will also assume the GPS will give us the Point of Aries location as
% well. (We never studied how this method works in depth, we will assume
% that it is less accurate and not as well defining.

% We will simulate noise by generating a little random vector that each
% sensor will mistake.

% Generate unit vectors for each sensor that they can be offset by.
% For GPS - use error = 0.017453292519943 radians maximum

% Errors
GPSError = 0.017453292519943 ; %Radians - Max for StarTracker
StarError = 3.878512000000000e-05 ; %Radians - Max for StarTracker

% Simulating Noise
ErrorVec = (rand(9,1)*2-1) ; %Randomized vector for all errors
MeasuredGPSError = ErrorVec(1:3,1)/norm(ErrorVec(1:3,1))*GPSError ; %Add to ideal meaasurement for actual (simulated)

MeasuredStarError1 = ErrorVec(4:6,1)/norm(ErrorVec(4:6,1))*StarError ; %Add to ideal for Startracker 1
MeasuredStarError2 = ErrorVec(7:9,1)/norm(ErrorVec(7:9,1))*StarError ; %Add to ideal for Startracker 2

%Use Errors to Generate Weights
WeightGPS = 1/(GPSError)^2 ; % Weight Factor for GPS
WeightStar = 1/(StarError)^2 ; % Weight Factor for StarTracker

% Ideal Measurements for each sensor (ASSUME THAT EACH SENSOR OUTPUTS BASED
% ON BODY ORIENTATION TO ECI)
% We will assume the sensor is pointed in the velocity direction of the
% orbit, other orietations are unknown, and thus assumed.
% Thus;
zBody = v/norm(v) ; %In ECI, need the x-axis in relation to this
xBody = cross(zBody,[4;5;2])/norm(cross(zBody,[4;5;2])) ; % Finds an Orthogonal Vector to the z-body, which is now  the x-body vector (Point of Aries Vector)
% GPS body relative to ECI
PointInECIGPS = [1;0;0] ; %Position of the point of Aries 
GPS = xBody ; %Ideal Position of Aries from Body

% Star Trackers body relative to ECI (Star Trackers will have same values
% (different error offsets)
PointInStarAries = [1;0;0] ; %Position of the point of Aries 
PointInStarZAxis = [0;0;1] ; %Z-axis of ECI
StarAries = xBody ; % Ideal Position of Aries from Body
StarZ = zBody ; %Uses the established z-orientation

% Measured Values Combine Ideals With Errors
MeasuredGPS = (GPS + MeasuredGPSError)/norm((GPS + MeasuredGPSError)) ;

MeasuredStarAries1 = (StarAries + MeasuredStarError1)/norm((StarAries + MeasuredStarError1)) ;
MeasuredStarZ1 = (StarZ + MeasuredStarError1)/norm((StarZ + MeasuredStarError1)) ;

MeasuredStarAries2 = (StarAries+MeasuredStarError2)/norm((StarAries+MeasuredStarError2)) ; 
MeasuredStarZ2 = (StarZ+MeasuredStarError2)/norm((StarZ+MeasuredStarError2)) ; 

% Now we must combine the sensor readings, the actual positions and the
% weights to get an initial orientation

% Known ECI
InputA = [PointInECIGPS,StarAries,StarZ,StarAries,StarZ] ;
% Sensor Readings
InputB = [MeasuredGPS,MeasuredStarAries1,MeasuredStarZ1,MeasuredStarAries2,MeasuredStarZ2] ;
% Weights
Weightings = [WeightGPS,WeightStar,WeightStar,WeightStar,WeightStar] ;
% Use Quest Function

[InitialOrientationECI] = Quest(InputA,InputB,Weightings) ;
% Use this initial Orientation for the detumble and continue from there.

%% Spacecraft Geometry  and MOI for the Detumble Phase
% ASSUME THAT THERE ARE NO PERTURBATIONS BECAUSE SATELLITE IS IN A CUBE,
% AND HOMOGENEOUS, ALL BUT MAGNETIC TORQUE ARE NEARLY 0 FOR A SATELLITE
% WITH EQUIVALENT MOI (ASSUME CUBE IS HOMOGENEOUS)

% Moments of Inertia
mass = 100 ; %kg
side = 2 ; %Length of sides in Meters
IAllCube = (mass*side^2)/6 ; %kg-m^2
ICube = eye(3).*IAllCube ; %kg-m^2

% Want to reach detumble in certain amount of time, but need thrusters that
% can provide the feedback. 
% Designing a control law for reaching detumble within 5 revolutions

% Control Law Constants-Detumble
zetaDet = .7 ; % Damping coefficient
SettlingTime = Period ; % Seconds - Settling Time


OmegaNat = 4.4/(SettlingTime*zetaDet) ; %radians/s -- Natural freq for the settling time and damping value

Kdeq =@(I,damp,OmegaNat) I*2*damp*OmegaNat ; % Function for finding Kd
Kpeq = @(I,OmegaNat) 2*I*OmegaNat^2 ; %Equation for each Kp Gain 
kpDetumble = Kpeq(ICube,OmegaNat) ; % Detumble Gain
kdDetumble = Kdeq(ICube,zetaDet,OmegaNat) ; % Detumble Gain

% Set Up Thrusters
% Total Of 12 Thrusters, 4 for each axis, counteracting and for both
% directions

% Each pair is .95 meters from the center of mass
ThrusterDistance = .95 ; %meters - This can be used to determine the thrust of each
% Know that 2 thrusters will be firing at once for each axis during control
% Assume their thrust can be increased and decreased in a threshold of
% design
ThrustEq = @(Torques) Torques./(2*ThrusterDistance) ; %Will give the thrust required from each individual thruster in each set
% Use this function to pull the control torque out of the statenew vector
Torquer = @(eps, omega , Kp , Kd) (-Kp*eps' - Kd*omega')' ; % Determine the control toque-N-m

%% Defining the Spacecraft Geometry and MOI For Operational Phase
magdi = [ 0 0 -0.5 ] ; % magnetic dipole moment A*m^2

%Calculating Moments of Inertia
mass = 100 ; %Total Mass in Kg
mw = 4.85 ;
panmass = .2*mass ; % mass of bus in kg
sensmass = .1*mass ; % sensor mass
busmass = mass - panmass - sensmass - 3*mw ; % solar panel mass
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
J = [ Ix 0 0 ; 0 Iy 0 ; 0 0 Iz ] ;

%Define Reaction Wheel Setup

rwheel = .247/2;
Lwheel = .085;
It = (1/2)*(mw*(3*rwheel^2+Lwheel^2)) ;
Is = (1/2)*(mw*rwheel^2) ;
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


% Control Law Constants
zeta = .7 ; % Damping coefficient
wn = .5 ;  % Wanted Natural Freq

kp = 2.*Isc.*wn^2 ;
kd = Isc.*2*zeta*wn ;

%Reaction Wheel Initial Angular Velocity
ww0 = [ 0 ; 0 ; 0 ] ; %Radians/s


%% Detumble Phase 
% Assume Omega =  [0.05; -0.02, -0.03] rad/s

% State Values
    % State 1 -> 4 are the quaternions, scalar first
    % State 5 -> 7 is the angular velocity of s/c
    % State 8 -> 10 is the position in ECI in km
    % State 11 -> 13 is the veloctiy in ECI in km/s
    % State 14 -> 16 is the Euler angles to ECI
    % State 17 -> 20 is the quaternions to LVLH
    % State 21 -> 23 is the Euler angles to LVLH
    % State 24 -> 27 is the quaternions for LVLH to ECI
    
%Use the sensor read values...
quatb2EciDet = dcm2quat(InitialOrientationECI)' ; %Finds initial body to ECI quaternion
%Initial Euler body to ECI
[Yaw0Det , Pitch0Det, Roll0Det]  = quat2angle(quatb2EciDet') ; %Finds initial Euler Angles
Anglesb2ECIDet = [Roll0Det , Pitch0Det, Yaw0Det]' ; %Places them into vector

% LVLH to ECI
Clvlhg = ECItoLVLH( r , v ) ;
quatLVLH2Eci = dcm2quat( Clvlhg )' ;

%Find Body to LVLH
quatb2LVLHDet = ComRotFind(quatb2EciDet',quatLVLH2Eci')' ;
[YawbL , PitchbL , RollbL ] = quat2angle(quatb2LVLHDet') ;
Anglesb2LVLHDet = [RollbL ; PitchbL ; YawbL] ;
% Set State
StateDet = [quatb2EciDet;OmegaDet;r;v;Anglesb2ECIDet;quatb2LVLHDet;Anglesb2LVLHDet;quatLVLH2Eci] ;
options = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
DetumbleDuration = SettlingTime*4 ;
[TimeDet,StateNewDet] = ode45(@DetumbleSolve,[0,DetumbleDuration],StateDet,options,ICube,muearth,kpDetumble,kdDetumble);

%Find the torques needed 
ThrustsNeeded = Thrusters(StateNewDet(:,5:7)',StateNewDet(:,1:4)',StateNewDet(:,8:10)',StateNewDet(:,11:13)',StateNewDet(:,24:27)',kpDetumble,kdDetumble,ThrusterDistance)' ; %Thrusts in Newtons

%% Operational Phase
%Find Initial Attitude relative to ECI (using end state of detumble)
% body to LVLH
eps = [ 0 0 0 ] ;
eta = 1 ;
quat0b = [ eta , eps ] ;
Cblvlh = quat2dcm( quat0b ) ;

% LVLH to ECI
Clvlhg = ECItoLVLH( StateNewDet(end,8:10)' , StateNewDet(end,11:13)' ) ;
quat0 = dcm2quat( Cblvlh*Clvlhg )' ;
%Initial Euler
[Yaw0 , Pitch0, Roll0]  = quat2angle(quat0') ; %Finds initial Euler Angles
Angles = [Roll0 , Pitch0, Yaw0]' ; %Places them into vector

% Establish State
    % State 1 -> 4 are the quaternions, scalar first
    % State 5 -> 7 is the angular velocity of s/c
    % State 8 -> 10 is the angular velocity of wheels
    % State 11 -> 13 is the position in ECI in km
    % State 14 -> 16 is the veloctiy in ECI in km/s
    % State 17 -> 19 is the Euler angles to ECI
    % State 20 -> 23 is the quaternions to LVLH
    % State 24 -> 26 is the Euler angles to LVLH
    % State 27 -> 30 is the quaternions for LVLH to ECI

State = [quat0;Omega;ww0;StateNewDet(end,8:10)';StateNewDet(end,11:13)';Angles;quat0b';[0;0;0];quat0] ;
options = odeset( 'RelTol' , 1e-8 , 'AbsTol' , 1e-8 ) ;
tend = 5*Period ;
torques = ones( 1 , 4 ) ;
% torques = [ 1 ; 1 ; 1 ; 0 ] ;
figure
hold on
[Time,StateNew] = ode45(@Operational,[0,tend],State,options,Isc,Iw,Areas,Positions,Norms,magdi,muearth,UTC,torques,0,kp,kd);
title( 'Magnitude of Torques on S/C with Upgrade Drag Model' )
xlabel( 'Time' )
ylabel( 'Torque (Nm)' )
legend( 'Drag Torque' , 'SRP' , 'Gravity Gradient' , 'Magnetic' )
hold off

%% Plots
% %% Plots for Detumble
% % Quaternions Body to LVLH
% TimeDet = TimeDet/(60*60) ;
% 
%     figure
%     hold on
%     plot(TimeDet,StateNewDet(:,17))
%     plot(TimeDet,StateNewDet(:,18))
%     plot(TimeDet,StateNewDet(:,19))
%     plot(TimeDet,StateNewDet(:,20))
%     title( 'Detumble: Quaternions Body to LVLH' )
%     legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' )
%     xlabel('Time [Hours]')
%     ylabel('Quaternion Value Magnitude')
%     hold off
%     
%     % Euler angles Body to LVLH
%     figure
%     hold on
%     plot(TimeDet,StateNewDet(:,21)*180/pi)
%     plot(TimeDet,StateNewDet(:,22)*180/pi)
%     plot(TimeDet,StateNewDet(:,23)*180/pi)
%     title( 'Detumble: Euler Angles Body to LVLH' )
%     legend( 'Roll' , 'Pitch' , 'Yaw' )    
%     xlabel('Time [Seconds]')
%     ylabel('Angle [Degrees]')
%     hold off
% 
% % Relative to ECI
%     % Quaternion LVLH to ECI
%     figure
%     hold on
%     plot(TimeDet,StateNewDet(:,24))
%     plot(TimeDet,StateNewDet(:,25))
%     plot(TimeDet,StateNewDet(:,26))
%     plot(TimeDet,StateNewDet(:,27))
%     title( 'Detumble: Quaternions LVLH to ECI' )
%     legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' , 'Orientation', 'Horizontal')
%     xlabel('Time [Hours]')
%     ylabel('Quaternion Value Magnitude')
%     hold off
%     
%     % Quaternion Body to ECI
%     figure
%     hold on
%     plot(TimeDet,StateNewDet(:,1))
%     plot(TimeDet,StateNewDet(:,2))
%     plot(TimeDet,StateNewDet(:,3))
%     plot(TimeDet,StateNewDet(:,4))
%     title( 'Detumble: Quaternions Body to ECI' )
%     legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' , 'Orientation', 'Horizontal')
%     xlabel('Time [Hours]')
%     ylabel('Quaternion Value Magnitude')
%     hold off
%     
%     % Euler angles Body to ECI
%     figure
%     hold on
%     plot(TimeDet,StateNewDet(:,14)*180/pi)
%     plot(TimeDet,StateNewDet(:,15)*180/pi)
%     plot(TimeDet,StateNewDet(:,16)*180/pi)
%     title( 'Detumble: Euler Angles Body to ECI' )
%     legend( 'Roll' , 'Pitch' , 'Yaw' , 'Orientation', 'Horizontal')    
%     xlabel('Time [Hours]')
%     ylabel('Angle [Degrees]')
%     hold off
%     
% 
% % Angular Velocity
%     figure
%     hold on
%     plot(TimeDet,StateNewDet(:,5))
%     plot(TimeDet,StateNewDet(:,6))
%     plot(TimeDet,StateNewDet(:,7))
%     title( 'Detumble: Angular Velocity of Spacecraft' )
%     legend( 'X-Axis' , 'Y-Axis' , 'Z-Axis'  )
%     xlabel('Time [Hours]')
%     ylabel('Angular Velocity [Radians/s]')
%     hold off
% 
% % Motor Thrusts Required
%     figure
%     for oo = 1:length(TimeDet)
%         NormThrust(oo) = norm(ThrustsNeeded(oo,:)) ;
%     end
%     CumulativeThrust = cumtrapz(TimeDet,NormThrust') ; %Sums the thrust required
%     hold on
%     plot(TimeDet,ThrustsNeeded(:,1)) ;
%     plot(TimeDet,ThrustsNeeded(:,2)) ;
%     plot(TimeDet,ThrustsNeeded(:,3)) ;
%     plot(TimeDet,NormThrust,'MarkerSize',2)
%     title('Thrust Required from Each Thruster (In a thruster quadfecta)')
%     xlabel('Time [Seconds]')
%     ylabel('Thrust [Newtons]')
%     legend( 'X-Axis Thruster' , 'Y-Axis Thruster' , 'Z-Axis Thruster' , 'Norm Thrust' )    
%     
%     figure
%     plot(TimeDet,CumulativeThrust,'MarkerSize',2)
%     title('Total Impulse Required From System of Thrusters ')
%     xlabel('Time [Hours]')
%     ylabel('Thrust [Newtons]')
% %     
% %% Plots For Operational 
% %Convert time to Hours
Time = Time/(60*60); %Seconds to Hours 
% % Relative to LVLH
% % Quaternion Body to LVLH
%     figure
%     hold on
%     plot(Time,StateNew(:,20))
%     plot(Time,StateNew(:,21))
%     plot(Time,StateNew(:,22))
%     plot(Time,StateNew(:,23))
%     title( 'Operational: Quaternions Body to LVLH' )
%     legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' )
%     xlabel('Time [Hours]')
%     ylabel('Quaternion Value Magnitude')
%     hold off
%     
%     % Euler angles Body to LVLH
%     figure
%     hold on
%     plot(Time,StateNew(:,24)*180/pi)
%     plot(Time,StateNew(:,25)*180/pi)
%     plot(Time,StateNew(:,26)*180/pi)
%     title( 'Operational: Euler Angle Body to LVLH' )
%     legend( 'Roll' , 'Pitch' , 'Yaw' )
%     xlabel('Time [Hours]')
%     ylabel('Angle [Degrees]')
%     hold off
%     upperlim = .001 ;
%     lowerlim = -upperlim ;
%     upper = upperlim*ones(length(Time),1) ;
%     lower = lowerlim*ones(length(Time),1) ;
%     figure
%     hold on
%     plot( Time , upper , 'k' , Time , lower , 'k' ) 
%     plot(Time,StateNew(:,24)*180/pi)
%     plot(Time,StateNew(:,25)*180/pi)
%     plot(Time,StateNew(:,26)*180/pi)
%     axis( [ 0 max(Time) -.002 .002 ] )
%     title( 'Operational: Euler Angle Body to LVLH' )
%     legend( 'Positive Pointing Error Limit' , 'Negative Pointing Error Limit' , 'Roll' , 'Pitch' , 'Yaw' )
%     xlabel('Time [Hours]')
%     ylabel('Angle [Degrees]')
%     hold off
% 
% % Relative to ECI
%     % LVLH to ECI
%     figure
%     hold on
%     plot(Time,StateNew(:,27))
%     plot(Time,StateNew(:,28))
%     plot(Time,StateNew(:,29))
%     plot(Time,StateNew(:,30))
%     title( 'Operational: Quaternions LVLH to ECI' )
%     legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)' , 'Orientation', 'Horizontal', 'Location','north')
%     xlabel('Time [Hours]')
%     ylabel('Quaternion Value Magnitude')
%     hold off
%     % Quaternion Body to ECI
%     figure
%     hold on
%     plot(Time,StateNew(:,1))
%     plot(Time,StateNew(:,2))
%     plot(Time,StateNew(:,3))
%     plot(Time,StateNew(:,4))
%     title( 'Operational: Quaternions Body to ECI' )
%     legend( 'Eta' , 'Eps(1)' , 'Eps(2)' , 'Eps(3)', 'Orientation', 'Horizontal', 'Location','north' )
%     xlabel('Time [Hours]')
%     ylabel('Quaternion Value Magnitude')
%     hold off
%     % Euler angles Body to ECI
%     figure
%     hold on
%     plot(Time,StateNew(:,17)*180/pi)
%     plot(Time,StateNew(:,18)*180/pi)
%     plot(Time,StateNew(:,19)*180/pi)
%     title( 'Operational: Euler Angle Body to ECI' )
%     legend( 'Roll' , 'Pitch' , 'Yaw' , 'Location','northwest')    
%     xlabel('Time [Hours]')
%     ylabel('Angle [Degrees]')
%     hold off
%     
% 
% % Angular Velocity
% figure
% hold on
% plot(Time,StateNew(:,5))
% plot(Time,StateNew(:,6))
% plot(Time,StateNew(:,7))
% title( 'Operational: Angular Velocity of Spacecraft' )
% legend( 'X-Axis' , 'Y-Axis' , 'Z-Axis' )
% xlabel('Time [Hours]')
% ylabel('Angular Velocity [Radians/s]')
% hold off
% 
% for ii = 1:length( Time )
%     wsc = StateNew( ii , 5:7 )' ;
%     h(1:3,ii) = Isc*wsc ;
% end
% 
% % Angular Momentum
% figure
% hold on
% plot(Time,h(1,:))
% plot(Time,h(2,:))
% plot(Time,h(3,:))
% title( 'Operational: Angular Momentum of Spacecraft' )
% legend( 'X-Axis' , 'Y-Axis' , 'Z-Axis' )
% xlabel('Time [Hours]')
% ylabel('Angular Momentum in Body[kg*m^2/s]')
% hold off
% 
% for ii = 1:length( Time )
%     Cbg = quat2dcm( StateNew( ii , 1:4 ) ) ;
%     heci(:,ii) = Cbg*h(:,ii) ;
% end
% % Angular Momentum
% figure
% hold on
% plot(Time,heci(1,:))
% plot(Time,heci(2,:))
% plot(Time,heci(3,:))
% title( 'Operational: Angular Momentum of Spacecraft' )
% legend( 'X-Axis' , 'Y-Axis' , 'Z-Axis' )
% xlabel('Time [Hours]')
% ylabel('Angular Momentum in ECI [kg*m^2/s]')
% hold off
% 
% Angular Velocity of Reaction Wheels
figure
hold on
plot(Time,StateNew(:,8)*(30/pi))
plot(Time,StateNew(:,9)*(30/pi))
plot(Time,StateNew(:,10)*(30/pi))
title( 'Operational: Angular Velocity of Reaction Wheels' )
legend( 'X-Axis Wheel' , 'Y-Axis Wheel' , 'Z-Axis Wheel' , 'Location','southeast' )
xlabel('Time [Hours]')
ylabel('Angular Velocity [RPM]')
hold off

% Angular Momentum of Reaction Wheels
for ii = 1:length( Time )
    w = StateNew(ii,8:10)' ;
    hw(:,ii) = Iw*w ;
end
figure
hold on
plot(Time,hw(1,:))
plot(Time,hw(2,:))
plot(Time,hw(3,:))
title( 'Operational: Angular Momentum of Reaction Wheels' )
legend( 'X-Axis Wheel' , 'Y-Axis Wheel' , 'Z-Axis Wheel' , 'Location','southeast' )
xlabel('Time [Hours]')
ylabel('Angular Momentum [kg*m^2/s]')
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


%% Quest Function
function [CbG] = Quest(InputVectrixA,InputVectrixB,Weight)
% input a vectrix of the input values, then in the same order, place the
% weight values. Outputs the rotation matrix
% Input Vectrix A is the known ECI or inertial vectors 
% Input Vectrix B are the found-by-device vectors

%First the B' matrix
B = zeros(3)' ; %Generates an empty matrix
for jj = 1:length(Weight)
    B = (B' + Weight(jj)*InputVectrixA(:,jj)*InputVectrixB(:,jj)')' ; %Defined in book as 25.7 
end

%Obtain S
S = B + B' ; %As seen in book
K12 = [(B(2,3)-B(3,2)),(B(3,1)-B(1,3)),(B(1,2)-B(2,1))]' ; %As defined on page 470
K22 = trace(B') ; %as defined in book at page 470

%Estimation of the lambda (eigenvalue is close to the sum of the weights)
Lam = sum(Weight) ; %Guess for lambda
e = inf ; %ensures step into while loop
if det(S) == 0
    a = K22^2 ; %Avoids the singularity of dividing by 0
else
    a = K22^2 - trace((det(S)*eye(3))/S) ; %defined for 20.20
end
b = K22^2 + K12'*K12 ;
c = det(S) + K12'*S*K12 ;
d = K12'*S^2*K12 ;
LamFun = @(L) L^4-(a+b)*L^2-c*L+(a*b+c*K22-d) ; %Function for Lambda
LamPrimeFun = @(L) 4*L^3 - (a+b)*2*L - c ; % Function for Lambda Prime
while e>1E-6
    LamPrev = Lam ;
    Lam = Lam - LamFun(Lam)/LamPrimeFun(Lam) ; %Newton Method Scheme for L
    e = abs(Lam-LamPrev) ; %Calculates the difference in between steps
end
if det(S) == 0
    Val = 0 ;
else 
    Val = trace((det(S)*eye(3))/S) ; 
end
alpha = Lam^2-K22^2 + Val ; %Alpha value
beta = Lam - K22 ; %Beta value
gamma = (Lam + K22)*alpha - det(S); %Gamma value
x = (alpha*eye(3)+beta*S+S.^2)*K12 ; %Xvector on 472
%Use 25.17  
q = 1/(sqrt(gamma^2+x'*x))*[x;gamma] ; %Finds the estimated quaternion from p (below 25.15)
CbG = quat2dcm([q(4);q(1:3)]') ; %converts estimate to a rotation matrix
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

%% Thruster Force Out For LVLH Control
function Thrust = Thrusters(wsc,quat,r,v,quatLVLH2eci,kpDetumble,kdDetumble,ThrusterDistance)
% Determine the total thrust of each thruster pair:
% Input all as 3 or 4 by x matrices
[~,M] = size(wsc) ;
Torquer = @(eps, omega , Kp , Kd) (-Kp*eps' - Kd*omega)' ; %Finds the torque for the control law
ThrustEq = @(Torques) Torques./(2*ThrusterDistance) ; %Will give the thrust required from each individual thruster in each set

for pp = 1:M
    Clvlhg = ECItoLVLH(r(:,pp),v(:,pp)) ;
    cquat = quatLVLH2eci(:,pp) ;
    wlvlh = Clvlhg*( cross( r(:,pp) , v(:,pp) ) / norm(r(:,pp))^2 ) ;
    werr = ( wsc(:,pp) - wlvlh ) ; 
    qerr = ComRotFind(quat(:,pp)',cquat') ;
    epse(1:3,1) = qerr( 2:4 ) ;
    TorquesForDetumble = Torquer(epse', werr , kpDetumble , kdDetumble) ; %Torques needed in N-m
    Thrust(:,pp) = ThrustEq(TorquesForDetumble) ; %Newtons - Thrust needed by each singular thruster (Two for Positive, Two for Negative) in each axis, Negative means opposite set of thrusters
end
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
r = r/norm(r) ;
sundirection = dot( r , s ) ;
theta = acosd(sundirection) ;
lam = asind( 6378 / norm(r) ) ; % https://ocw.tudelft.nl/wp-content/uploads/AE2104-Orbital-Mechanics-Slides_10.pdf
eclipseang = 180-lam/2 ;
if theta <= eclipseang 
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

%% General ODE Functions
%% DETUMBLE
function dstate = DetumbleSolve(t,state,Isc,mu,kp,kd)
% State 1 -> 4 are the quaternions, scalar first
% State 5 -> 7 is the angular velocity of s/c
% State 8 -> 10 is the position in ECI in km
% State 11 -> 13 is the veloctiy in ECI in km/s
% State 14 -> 16 is the Euler angles to ECI
% State 17 -> 20 is the quaternions to LVLH
% State 21 -> 23 is the Euler angles to LVLH
% State 24 -> 27 is the quaternions for LVLH to ECI

% Remove before turning in, just for seeing progress
tdetumble = t

    dstate = zeros(26,1) ;
% Set up names
    Ix = Isc(1,1) ;
    Iy = Isc(2,2) ;
    Iz = Isc(3,3) ;

    quat = state( 1:4 )' ;
    wsc = state( 5:7 ) ;
    v = state( 11:13 ) ;
    r = state( 8:10 ) ;
    quatLVLH = state(17:20)' ;
    quatLVLH2eci = state(24:27)' ;
    
% Define Rotation based off of quat
    eta = quat( 1 ) ;
    eps = quat( 2:4 )' ;
    etaLVLH = quatLVLH( 1 ) ;
    epsLVLH = quatLVLH( 2:4 )' ;
    etaLVLH2eci = quatLVLH2eci( 1 ) ;
    epsLVLH2eci = quatLVLH2eci( 2:4 )' ;    
    
    if t >= 20000
        tl = t
    end
    
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
    dstate(17:20) = [ etadote , epsdote' ]' ;   
% change quaternions LVLH to ECI
    epsdotl2e = .5 * ( etaLVLH2eci*eye(3)*wlvlh + cross( epsLVLH2eci , wlvlh ) ) ;
    etadotl2e = -.5 * epsLVLH2eci' * wlvlh ;
    dstate(24:27) = [ etadotl2e , epsdotl2e' ]' ;   

%Simplified Orbital Solver, dState 11 -> 13 is the acceleration 
    dstate(11:13) = -mu*r/norm(r)^3 ; %Position to Acceleration Function
    dstate(8:10) = v ; %velocity in velocity out

%Define Control Torque
    cTorque = -kp*epse - kd*werr ;
    
%Use Torques to find the proper omegadot in body frame
    wcross = [ 0 -wsc(3) wsc(2) ; wsc(3) 0 -wsc(1) ; -wsc(2) wsc(1) 0 ] ;
    dwsc = (Isc)\( cTorque - wcross*Isc*wsc ) ;
    dstate(5:7) = dwsc ; %Inputs the omega dot due to given torques and motion

% Derivitive for ECI Euler
    roll = (state(14)) ;
    pitch = (state(15)) ;
    dstate(14:16) = [1,sin(roll)*tan(pitch),cos(roll)*tan(pitch);0,cos(roll),-sin(roll);0,sin(roll)*sec(pitch),cos(roll)*sec(pitch)]*wsc ; %Euler angle change rate ; % Euler Body to ECI

% Derivative for Body to LVLH Euler
    rollLVLH = (state(21)) ;
    pitchLVLH = (state(22)) ;
    dstate(21:23) = [1,sin(rollLVLH)*tan(pitchLVLH),cos(rollLVLH)*tan(pitchLVLH);0,cos(rollLVLH),-sin(rollLVLH);0,sin(rollLVLH)*sec(pitchLVLH),cos(rollLVLH)*sec(pitchLVLH)]*werr ; %uiler angle change rate ; % Euler Body LVLH

end
%% OPERATIONAL PHASE
function dstate = Operational(t,state,Isc,Iw,Areas,Positions,Norms,MagDi,mu,UTC,Torques,Upgrade,kp,kd)
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

% Derivative for Body to LVLH Euler
rollLVLH = (state(24)) ;
pitchLVLH = (state(25)) ;
dstate(24:26) = [1,sin(rollLVLH)*tan(pitchLVLH),cos(rollLVLH)*tan(pitchLVLH);0,cos(rollLVLH),-sin(rollLVLH);0,sin(rollLVLH)*sec(pitchLVLH),cos(rollLVLH)*sec(pitchLVLH)]*werr ; %uiler angle change rate ; % Euler Body LVLH

end



