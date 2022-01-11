%Pre Lab 7
%Liam Hood

% y'' + y*y' + y = 0
% y(1)' = y(2) = y'
% y(2)' = -y(1) - y(1)*y(2)

%Evaluates the system over the tspan seconds with the initial values given
%by y0. The first row of y0 is the initial y value while the second row is
%the initial value of y'
tspan = [ 0 10 ] ;
y0 = [ 1 ; 0 ] ;
[ t , y ] = ode45( @part4 , tspan , y0 ) ;

figure
plot( t , y )
title( 'IVP with ODE45' )
xlabel( 'Time' )
ylabel( 'y values' )
legend( 'y' , 'y''' ) 
grid

function dydt = part4( t , y )
%Sets up the system of equations in a form that ode45 can use. The first is
%y' and the second row is y''
    dydt =  [ y(2) ; -y(1) - y(1)*y(2) ] ;
end

