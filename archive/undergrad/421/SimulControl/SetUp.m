clear ; close all ; clc ;
% COES and initial state
h = 53335.2 ; % km^2/s
ecc = 0 ; 
RAAN = 0 ;
inc = 98.43 ; % degrees
omega = 0 ;
theta = 0 ;
mu = 398600 ;
[ r0 , v0 ] = coes2state(  h , ecc , theta , RAAN , omega , inc , mu ) ;

% Areas
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

%Calculating Moments of Inertia
side = 2 ; % side length in m
busmass = 500 ; % mass of bus in kg
magdi = [ 0 0 -0.5 ] ; % magnetic dipole moment A*m^2
sensmass = 100 ; % sensor mass
panmass = 20 ; % solar panel mass
% Geometry means that x and y center is aligned with the center of mass.
% There is an offset for the z axis, as the sensor pulls it off the
% Z offset calc
zoff = sensmass*1.5/(sensmass+2*panmass+busmass) ;
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
I = [ Ix 0 0 ; 0 Iy 0 ; 0 0 Iz ] ;

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
Norms = [nx     nx      -nx     -nx     ny      -ny     nz      nz      -nz     -nz     nx      -nx     ny      -ny     nz      nx      -nx     ny      -ny     nz      -nz     ] ; %Normal Vectors for each face
% The full positions and geometries of the faces are defined

%Omega Vector Initial
w0 = [0;-0.001047;0] ; %rad/s - This is in body frame

% Sun Vector
s = [ -1 0 0 ] ;

% Start time 
UTC = [ 2019 , 3 , 19 , 12 , 0 , 0 ] ;

%Find Initial Attitude relative to ECI
% body to LVLH
eps0b = [ 0 0 0 ] ;
eta0b = 1 ;
quat0b = [ eta0b , eps0b ] ;
Cblvlh = quat2dcm( quat0b ) ;

% LVLH to ECI
Clvlhg = ECItoLVLH( r0 , v0 ) ;
quat0 = dcm2quat( Cblvlh*Clvlhg )' ;
eta0 = quat0(1) ;
eps0 = quat0(2:4) ;
%Initial Euler
[Yaw0 , Pitch0, Roll0]  = quat2angle(quat0') ; %Finds initial Euler Angles
Angles = [Yaw0 , Pitch0, Roll0]' ; %Places them into vector






function [ r , v ] = coes2state(  h , ecc , theta , RAAN , omega , inc , mu )
    r_peri = (h^2/mu) * ( 1/( 1 + ecc*cosd(theta) ) ) * [ cosd( theta ) ; sind( theta ) ; 0 ] ;
    v_peri = (mu/h) * [ -sind( theta ) ; ecc+cosd(theta) ; 0 ] ;

    d2r = pi/180 ;
    RAAN = d2r*RAAN ;
    omega = d2r*omega ;
    inc = d2r*inc ;
    Q(1,1) = -sin(RAAN)*cos(inc)*sin(omega) + cos(RAAN)*cos(omega) ;
    Q(1,2) = -sin(RAAN)*cos(inc)*cos(omega) - cos(RAAN)*sin(omega) ;
    Q(1,3) = sin(RAAN)*sin(inc) ;
    Q(2,1) = cos(RAAN)*cos(inc)*sin(omega) + sin(RAAN)*cos(omega) ;
    Q(2,2) = cos(RAAN)*cos(inc)*cos(omega) - sin(RAAN)*sin(omega) ;
    Q(2,3) = -cos(RAAN)*sin(inc) ;
    Q(3,1) = sin(inc)*sin(omega) ;
    Q(3,2) = sin(inc)*cos(omega) ;
    Q(3,3) = cos(inc) ;

    r = Q*r_peri ;
    v = Q*v_peri ;
end

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