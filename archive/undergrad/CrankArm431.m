clear ; close all ; clc ;
n = 1e3 ; % number of steps
L = .16226 ; % total length in meters
F = 100*9.81*6 ; % Force 220kg person at 6g
x = linspace( 0 , L , n ) ; % create x vector of n steps over total length
bmin = .015 ; % base of rectangle
bmax = .015 ; 
hmin = .014 ; % minimum height
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


mshearnt = 3*F/( max(b)*max(h) + pi*(max(r)^2) ) ;
disp([ 'The maximum shear stress is along the sides in the center and is  ', num2str( mshearnt*1e-6 ) , ' MPa, along the whole length in a straight crank' ])
shearnt = mshearnt*ones(1,n) ;
shear = 3*F./( b.*h + pi.*(r.^2) ) ;
[ mshear , index ] = max( shear ) ;
xshear = x( index ) ;
disp([ 'The maximum shear stress is along the sides in the center and is ', num2str( mshear*1e-6 ) , ' MPa, ' , num2str( xstress ) , ' meters from pedal, in a tapered crank' ])
figure 
plot( x , stressnt*1e-6 , x , stress*1e-6 , x , shearnt*1e-6 , x , shear*1e-6 )
title( 'Stress on a Beam' )
xlabel( 'Distance from Pedal (m)' )
ylabel( 'Stress (MPa)' )
legend( 'Normal Stress, Straight Arm' , 'Normal Stress, Tapered Arm' , 'Shear Stress, Straight Arm' , 'Shear Stress, Tapered Arm' , 'Location' , 'east' )
    
% Factor of Safety
FS = 1.1 ;

% Taper
    % Tresca 
    shearmax = mstress / 2 ;
    shearyt = FS*shearmax ;
    criticalT = shearyt*2 ;
    disp([ 'Critical Stress by Tresca for tapered is ' , num2str( criticalT*1e-6 ) , ' MPa for a Factor of safety of ' , num2str(FS) ]) ;
    
    % Von Meises
    qmax = sqrt( 2*mstress^2 ) ;
    qcrit = FS*qmax ;
    criticalVM = sqrt( ( qcrit^2 )/2 ) ;
    disp([ 'Critical Stress by Von Meises for tapered is ' , num2str( criticalVM*1e-6 ) , ' MPa for a Factor of safety of ' , num2str(FS) ]) ;
    
% No Taper  
    % Tresca 
    shearmaxnt = mstressnt / 2 ;
    shearytnt = FS*shearmaxnt ;
    criticalTnt = shearytnt*2 ;
    disp([ 'Critical Stress by Tresca for non-tapered is ' , num2str( criticalTnt*1e-6 ) , ' MPa for a Factor of safety of ' , num2str(FS) ]) ;
    
    % Von Meises
    qmaxnt = sqrt( 2*mstressnt^2 ) ;
    qcritnt = FS*qmaxnt ;
    criticalVMnt = sqrt( ( qcritnt^2 )/2 ) ;
    disp([ 'Critical Stress by Von-Meises for non-tapered is ' , num2str( criticalVMnt*1e-6 ) , ' MPa for a Factor of safety of ' , num2str(FS) ]) ;