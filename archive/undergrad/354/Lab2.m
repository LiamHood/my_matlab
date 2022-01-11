clear; close all; clc;
load( 'Langmuir_practice.txt' )
V = Langmuir_practice(:,1) ;
I = Langmuir_practice(:,2)*(10^-3) ;
[ ioffset , index ] = min(I) ;
I = I - ioffset ;
Vfloat = V( index ) ;
    eps0 = 8.85418782e-12 ;

e = 1.60217662e-19 ;
me = 9.1938356e-31;
mi = 6.65e-26 ;
k = 1.38064852e-23 ;
Te = e*(V(500)-V(390))/(k*log(I(500)/I(390))) ; % kelvin
Ti = Te ;
% Aprobe = I./j ;
Aprobe = 1 ;
veth = sqrt( 8*k*Te/me )*Aprobe;
vith = sqrt( k*Te/mi )*Aprobe;
n = -ioffset./(.25*e*vith*Aprobe) ;
j = n*e*veth ;
    debye = sqrt( eps0*k*Te/(n*(e^2)) ) ;
    pp = (4/3)*pi*debye^3*n ;
freq = 8.98*sqrt(n) ;
slope = (I(500)/I(390))/(V(500)-V(390)) ;
midline = slope.*V - I(400) ;
Ies = (1/4)*e*n*veth*Aprobe ;    
PlasPot = ( Ies + slope*V(445) - I(445) )/slope ;
figure
semilogy( V , I )



