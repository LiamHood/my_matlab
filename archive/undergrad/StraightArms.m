clear ; close all ; clc ;
n = 1e3 ; % number of steps
L = .16226 ; % total length in meters
F = 100*9.81*6 ; % Force 220kg person at 6g
x = linspace( 0 , L , n ) ; % create x vector of n steps over total length
bmin = .015 ; % base of rectangle
bmax = .015 ; 
hmin = .040 ; % minimum height
hmax = .040 ; % maximum height
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


stresslarge = (M.*c)./I ; % Stress of tapered beam
[ mstresslarge , index ] = max( stresslarge ) ;
xstresslarge = x( index ) ;
disp([ 'The maximum stress is ', num2str( mstresslarge*1e-6 ) , ' MPa, ' , num2str( xstresslarge ) , ' meters from pedal, in a tapered crank' ])

    bmin = .015; % base of rectangle
    bmax = .015 ; 
    hmin = .035 ; % minimum height
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
    
stress = (M.*c)./I ; % Stress of tapered beam
[ mstress , index ] = max( stress ) ;
xstress = x( index ) ;
disp([ 'The maximum stress is ', num2str( mstress*1e-6 ) , ' MPa, ' , num2str( xstress ) , ' meters from pedal, in a tapered crank' ])

    bmin = .015; % base of rectangle
    bmax = .015 ; 
    hmin = .030 ; % minimum height
    hmax = .030 ; % maximum height
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
    
stresssmall = (M.*c)./I ; % Stress of tapered beam
[ mstresssmall , index ] = max( stresssmall ) ;
xstresssmall = x( index ) ;
disp([ 'The maximum stress is ', num2str( mstresssmall*1e-6 ) , ' MPa, ' , num2str( xstresssmall ) , ' meters from pedal, in a tapered crank' ])

figure 
plot( x , stresslarge*1e-6 , x , stress*1e-6 , x , stresssmall*1e-6 )
title( 'Normal Stress in Straight Crank Arms of Varying Heights' )
xlabel( 'Distance from Pedal (m)' )
ylabel( 'Stress (MPa)' )
legend( '40mm' , '35mm' , '40mm' ,  'Location' , 'northwest' )
    
