clear ; close all; clc;
disp( 'For a rectangular prism s/c' )
Qdis = 200 ; % Heat dissipated
rad = 1 ; % radius in AU
A = 2*2*6 ; % total surface area in square meters
Awet = 1/6 * A ; % lit surface area in square meters
S = 1370 ; % Energy from sun in W/m^2
sbc = 5.67e-8 ; % stefan boltzmann constant
alpha = linspace(0,1,1e2) ; % range of possible absorbance
Tc = -5 ; % Temperature in celsius
T = 273.15 + Tc ; % Temp in kelvin
Tch = 50 ; % Temperature in celsius
Th = 273.15 + Tch ; % Temp in kelvin

Qsun = ( S / ( rad^2 )) * Awet * alpha ; % Heat from sun
intermediate = ( Qsun + Qdis )/(sbc*T^4) ;
eps = intermediate/A ; % emmissivity as a function of alpha
intermediateh = ( Qsun + Qdis )/(sbc*Th^4) ;
epsh = intermediateh/A ; % emmissivity as a function of alpha
figure
plot( alpha , eps , alpha , epsh )
xlabel( 'Absorbtance' )
ylabel( 'Emmissivity' )
title( 'Surface properties of coating' )
legend( 'To maintain cold temperature' , 'To maintain hot temperature' )

disp( 'bare aluminum will work as a coating at its highest alpha and epsilon' )
disp( 'to keep the s/c within operating temp at 1 AU' )
disp( ' ' )

alpha2 = .17 ; % using aluminumm
eps2 = .1 ; % using aluminum
rad2 = 19.2 ; % radius in AU
Tc2 =  5 ; % temp to maintain in Celsius
T2 = 273.15 + Tc2 ; % temp in kelvin

Qsun2 = ( S / ( rad2^2 )) * Awet * alpha2 ; % heat from sun at Uranus
Qout2 = A*eps2*sbc*T2^4 ; % Heat leaving due to radiation
Qheat2 = Qout2 - Qsun2 - Qdis ; % heat heater needs to generate in watts
disp([ 'Need a heater of more than ' ,  num2str(Qheat2) , ' watts to keep above 5 degrees C at Uranus' ])
disp( ' ' )

%% Sphere
disp( 'For a spherical s/c' )
A2 = 4*pi ;
Awet2 = 1/2 * A2 ;

Qsun3 = ( S / ( rad^2 )) * Awet2 * alpha ; % heat from sun
intermediate3 = ( Qsun3 + Qdis )/(sbc*T2^4) ;
eps3 = intermediate3/A2 ; % emissitivity
intermediate3h = ( Qsun3 + Qdis )/(sbc*Th^4) ;
eps3h = intermediate3h/A2 ;
figure
plot( alpha , eps3 , alpha , eps3h )
xlabel( 'Absorbtance' )
ylabel( 'Emmissivity' )
title( 'Surface properties of coating on Spherical S/C' )
legend( 'To maintain cold temperature' , 'To maintain hot temperature' , 'Location' , 'northwest' )

disp( 'Any aluminized kapton will keep the s/c within temp range at Earth' )
disp( ' ' )

alpha4 = .38 ;
eps4 = .67 ;

Qsun4 = ( S / ( rad2^2 )) * Awet2 * alpha4 ;
Qout4 = A2*eps4*sbc*T2^4 ;
Qheat4 = Qout4 - Qsun4 - Qdis ;

disp([ 'Need a heater of more than ' ,  num2str(Qheat4) , ' watts to keep above 5 degrees C at Uranus' ])


