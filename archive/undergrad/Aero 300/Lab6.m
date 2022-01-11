%Liam Hood
%Aero 300
%Lab 6

clear
close all

%% 1. Numeric Differentiation
displacement = load( 'dispData.mat' );
n = length( displacement.y ) ;
h = displacement.t(2) ;

% Two-point
v2p = zeros( 1 , 100 ) ;
for ii = 1:n-1 %Runs for all points save for the last data point 
    v2p(ii) = ( displacement.y( ii+1 ) - displacement.y( ii ) ) / h ; %Formula for two-point differentiation
end

% Three-point
v3p = zeros( 1 , 100 ) ;
for ii = 2:n-1
    v3p(ii) = ( displacement.y( ii+1 ) - displacement.y( ii-1 ) ) / ( 2*h ) ; %Formula for three-point differentiation
end

% Three point acceleration
a3p = zeros( 1 , 100 ) ;
for ii = 2:n-1
    a3p(ii) = ( displacement.y( ii-1 ) - 2*displacement.y( ii ) +  displacement.y( ii+1 ) ) / ( h^2 ) ; %Formula for three-point second derivative
end

% Dacksolved Displacement by kinematics
dExp = [ 0 , zeros( 1 , 99 ) ] ;
for ii = 2:n
    dExp(ii) = dExp(ii-1) + v3p(ii-1)*h ;
end

%Best fit for displacement
dPc = polyfit( displacement.t , displacement.y , 2 ) ;
dP = dPc(1)*( displacement.t.^2 ) + dPc(2).*displacement.t + dPc(3) ;

%Derivative of displacement
vPc = polyder( dPc ) ;
vP = vPc(1)*displacement.t + vPc(2) ;

%derivative of velocity
aPc = polyder( vPc ) ;
aP = linspace( aPc , aPc , n ) ;

% Plots
figure( 'OuterPosition' , [ 0 0 1200 800 ] ) 

subplot( 2 , 3 , 1 ) %Top left
plot( displacement.t , displacement.y , '.' ) % measured displacement
title( 'Measured Displacement vs Time' ) % title of plot
xlabel( 'Time (s)' ) % title of x axis
ylabel( 'Displacement(m)' ) % title of y axis
grid %displays grid on graph

subplot( 2 , 3 , 4 ) %Bottom left
plot( displacement.t , dExp , '.' )  %Displacement by best fit
title( 'Expected Displacement vs Time' )
xlabel( 'Time (s)' )
ylabel( 'Displacement (m)' )
grid

subplot( 2 , 3 , 2 ) %top center
plot( displacement.t , v2p , '.' ) %Velocity from two-point differentiation
title( 'Velocity (2-point) vs Time' )
xlabel( 'Time (s)' )
ylabel( 'Velocity (m/s)' )
grid

subplot( 2 , 3 , 5 ) %bottom center
plot( displacement.t , v3p , '.' ) %Velocity from three-point differentiation
title( 'Velocity (3-point) vs Time' )
xlabel( 'Time (s)' )
ylabel( 'Velocity (m/s)' )
grid

subplot( 2 , 3 , 3 ) %top right
plot( displacement.t , a3p , '.' ) %Acceleration from three-point double differentiation
title( 'Acceleration vs Time' )
xlabel( 'Time (s)' )
ylabel( 'Acceleration (m/s^s)' )
grid

 %%  Better Graphs for 1
 figure( 'OuterPosition' , [ 0 0 800 1000 ] ) 
 subplot( 3 , 1 , 1 ) %Top  
 plot( displacement.t , displacement.y , 'b.-' , displacement.t , dExp , 'r.-' , displacement.t , dP , 'k-' ) %displacement
 title( 'Displacement vs Time' ) % title of plot
 xlabel( 'Time (s)' ) % title of x axis
 ylabel( 'Displacement(m)' ) % title of y axis
 grid %displays grid on graph
 legend( 'Measured Displacement' , 'Expected Displacement Velocity and Kinematics' , 'Expected Displacement Based on Best Fit Curve' )
 
 

 subplot( 3 , 1 , 2 ) %center
 plot( displacement.t , v2p , 'b.-' , displacement.t , v3p , 'r.-' , displacement.t , vP ) %Velocity
 title( 'Velocity vs Time' )
 xlabel( 'Time (s)' )
 ylabel( 'Velocity (m/s)' )
 grid
 legend( '2-point' , '3-point' , 'From kinematics and best fit' )
 
 
 subplot( 3 , 1 , 3 ) %bottom
 plot( displacement.t , a3p , '.-' , displacement.t , aP ) %Acceleration
 title( 'Acceleration vs Time' )
 xlabel( 'Time (s)' )
 ylabel( 'Acceleration (m/s^s)' )
 legend( 'Acceleration by three-point method' , 'Acceleration by best fit and kinematics' )
 grid

disp( 'Taking the numerical derivative of noisy data leads to even noiser ' )
disp( 'derivative estimates. A better way to find the acceleration from this' )
disp( 'data would be to fit a curve to the displacement data and take the ' )
disp( 'second derivative of that.' )

%% 2. Integration
accel = load( 'accelData.mat' );
n = length( accel.acc );
h = accel.t(2) ;

% Left-Handed Rectangular
vR = zeros( 1 , n ) ;
for ii = 2:n
    vR(ii) = vR(ii-1) + h*accel.acc(ii-1) ; %Left hand rectangular formula for velocity
end

dR = zeros( 1 , n ) ;
for ii = 2:n
    dR(ii) = dR(ii-1) + h*vR(ii-1) ; %Uses velocity and left hand rectangular to find displacement
end

%Trapezoid 
vT = zeros( 1 , n ) ;
for ii = 2:n
    vT(ii) = vT(ii-1) + (h/2)*( accel.acc(ii-1) + accel.acc(ii) ) ; %Trapezoid method for velocity
end

dT = zeros( 1 , n ) ;
for ii = 2:n
    dT(ii) = dT(ii-1) + (h/2)*( vT(ii-1) + vT(ii) ) ; %Velocity and trapezoid method to find displacement
end

%Simpson's 
vS = zeros( 1 , n ) ;
for ii = 2:n-1
    vS(ii) = vS(ii-1) + (h/6)*( accel.acc(ii-1) + 4*accel.acc(ii) + accel.acc(ii+1) ) ; %Simpson's method to find velocity
end

dS = zeros( 1 , n ) ;
for ii = 2:n-1
    dS(ii) = dS(ii-1) + (h/6)*( vS(ii-1) + 4*vS(ii) + vS(ii+1) ) ; %Simpsons's and velocity to find displacement
end

% Expected
aEc = polyfit( accel.t' , accel.acc , 0 ); %Best fit acceleration
aE = linspace( aEc , aEc , n ) ;

vE = aEc*accel.t ; %Velocity predicted by kinematics

dE = .5 * aEc *( accel.t .^ 2 ) ; %Displacement predicted by kinematics


figure( 'OuterPosition' , [ 0 0 600 1000 ] )
subplot( 3 , 1 , 1 )
plot( accel.t , accel.acc , '.-b' , accel.t , aE , '-r' ) %Acceleration plot
title( 'Acceleration vs Time' )
ylabel( 'Acceleration (m/s^2)' )
xlabel( 'Time (s)' )
legend( 'Measured Acceleration' , 'Acceleration Best Fit' )

subplot( 3 , 1 , 2 )
plot( accel.t , vR , '.-r' , accel.t , vT , '.-b' , accel.t , vS , '.-k' , accel.t , vE , '-g') %Velocity plot
title( 'Velocity vs Time' )
ylabel( 'Velocity (m/s)' )
xlabel( 'Time (s)' )
legend( 'Left-Hand Rectangular' , 'Trapezoid Method' , 'Simpson''s' , 'Kinematics' )

subplot( 3 , 1 , 3 )
plot( accel.t , dR , '.-r' , accel.t , dT , '.-b' , accel.t , dS , '.-k' , accel.t , dE , '-g' ) %Acceleration plot
title( 'Displacement vs Time' )
ylabel( 'Displacement (m)' )
xlabel( 'Time (s)' )
legend( 'Left-Hand Rectangular' , 'Trapezoid Method' , 'Simpson''s' , 'Kinematics' )

%% 3. Discussion
