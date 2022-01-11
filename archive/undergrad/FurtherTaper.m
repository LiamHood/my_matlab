clear ; close all ; clc ;
n = 1e3 ; % number of steps
L = .16226 ; % total length in meters
F = 100*9.81*6 ; % Force 220kg person at 6g
x = linspace( 0 , L , n ) ; % create x vector of n steps over total length
bmin = .015 ; % base of rectangle
bmax = .015 ; 
hmin = .0155 ; % minimum height
hmax = .035 ; % maximum height
b = (( bmax - bmin )/L).*x + bmin ; % base function
r = b/2 ; % radius of half circle
h = (( hmax - hmin )/L).*x + hmin ; % height function
M = F*x ; % Moment function
c = (h./2) + r ; % c function
for ii = 1:n
Icirc(ii) = 2*( pi/8 - 8/(9*pi) )*r(ii)^4 ; % Inertia of half circle
Irec(ii) = (1/12)*b(ii)*h(ii)^3 ; % inertia of rectangle
Ioffset(ii) = r(ii)^2*pi*( (4*r(ii))/(3*pi) + (h(ii)/2) )^2 ; % inertia to offset circle
I(ii) = Irec(ii) + Icirc(ii) + Ioffset(ii) ; % total inertia function
end
stressnt = F.*x*max(c)/max(I) ; % Stress of non tapered beam
[ mstressnt , indexnt ] = max( stressnt ) ;
xstressnt = x( indexnt ) ;
disp([ 'The maximum stress is ', num2str( mstressnt*1e-6 ) , ' MPa, ' , num2str( xstressnt ) , ' meters from pedal, in a non-tapered crank' ])


stress = (M.*c)./I ; % Stress of tapered beam
[ mstress , index ] = max( stress ) ;
xstress = x( index ) ;
disp([ 'The maximum stress is ', num2str( mstress*1e-6 ) , ' MPa, ' , num2str( xstress ) , ' meters from pedal, in a tapered crank' ])

    bmin = .00665 ; % base of rectangle
    bmax = .015 ; 
    hmin = .0155 ; % minimum height
    hmax = .035 ; % maximum height
    b = (( bmax - bmin )/L).*x + bmin ; % base function
    r = b/2 ; % radius of half circle
    h = (( hmax - hmin )/L).*x + hmin ; % height function
    M = F*x ; % Moment function
    c = (h./2) + r ; % c function
    for ii = 1:n
    Icirc(ii) = 2*( pi/8 - 8/(9*pi) )*r(ii)^4 ; % Inertia of half circle
    Irec(ii) = (1/12)*b(ii)*h(ii)^3 ; % inertia of rectangle
    Ioffset(ii) = r(ii)^2*pi*( (4*r(ii))/(3*pi) + (h(ii)/2) )^2 ; % inertia to offset circle
    I(ii) = Irec(ii) + Icirc(ii) + Ioffset(ii) ; % total inertia function
    end
    
stressft = (M.*c)./I ; % Stress of further tapered beam
[ mstressft , index ] = max( stressft ) ;
xstressft = x( index ) ;
disp([ 'The maximum stress is ', num2str( mstressft*1e-6 ) , ' MPa, ' , num2str( xstressft ) , ' meters from pedal, in a tapered crank' ])
figure 
plot( x , stressnt*1e-6 , x , stress*1e-6 , x , stressft*1e-6 )
title( 'Normal Stress' )
xlabel( 'Distance from Pedal (m)' )
ylabel( 'Stress (MPa)' )
legend( 'Straight Crank Arm' , 'Mild Taper' , 'Extreme Taper' ,  'Location' , 'southeast' )
 
