% Lab 7
% Liam Hood

close all ;
clear ;

%% 1: Simple
%Givens
fun = @(t,y) ( y - t - 1 )^2 + 2 ;
tspan = [ 0 pi/3 ] ;
y0 = 1 ;
h = .01 ;
rTol = 10^-6 ;

%Solve ODE
[t,y] = rkf45(fun, tspan, y0, h, rTol) ; %Run my function

[ tODE , yODE ] = ode45( fun , tspan , y0 ); %Run built in function

yAct = @(t) tan(t) + t + 1 ; %Analytical solution

figure %Plot of ODE solvers and actual solution
plot( t , y , '.b' , tODE , yODE , '.r' , tODE , yAct( tODE ) , 'k' ) 
ylabel( 'y' )
xlabel( 'x' )
title( 'RKF45' )
legend( 'My RKF45' , 'MATLAB ODE45' , 'Actual Solution' )

% Error
errorRKF = abs( y - yAct( t ) ) ; %Error in my function
errorODE = abs( yODE - yAct( tODE ) ) ; %Error in Matlab function

figure
%RKF Error
subplot( 2 , 1 , 1 )
semilogy( t , errorRKF , '.-b' ) %Y axis is logarithmic scale
title( 'Error of My RKF Compared to Function' )
xlabel( 'x' )
ylabel( 'Log of Error' )
%ODE Error
subplot( 2 , 1 , 2 )
semilogy( tODE , errorODE , '.-r' ) %Y axis is logarithmic scale
title( 'Error of Matlab ODE45 Compared to Function' )
xlabel( 'x' )
ylabel( 'Log of Error' )

%% 2: Lorenz
%Inputs
fun = @lorenz ;
interval = [ 0 50 ] ;
y0 = [ -20 , 30 , 30 ] ;
h = .01 ;
rTol = 10^-1 ;

%Evaluate lorenz
tic
[ t , pos ] = rkf45( fun , interval , y0 , h , rTol ) ; %RKF with wide tolerance
t1 = toc ;
tic
[ tM , posM ] = ode45( fun , interval , y0 ); %ODE using default tolerances
t2 = toc ;
tic
[ t_3 , pos3 ] = rkf45( fun , interval , y0 , h , 10^-15 ) ; %Tight tolerance RKF
t3 = toc ;

%Plot of different lorenz solution
figure( 'pos' , [ 200 , 100 , 900 , 600 ] )
plot3( pos(1,:) , pos(2,:) , pos(3,:) , posM(:,1) , posM(:,2) , posM(:,3) , pos3(1,:) , pos3(2,:) , pos3(3,:))
legend( 'Tolerance of .1' , 'ODE45 Solution' , 'Tolerance of 10^-10' )
title( 'Lorenz' )
xlabel( 'x' ) 
ylabel( 'y' )
zlabel( 'z' )

disp( 'The 3d plots don''t show up when I publish' )

%% Functions

function [tt,y] = rkf45(fun, tspan, y0, h, rTol)
%Runge-Katta_Fehlberg order 4/5. Input function handle fun; time span as a 
%vector tspan [ lower upper ] ; initial y value y0; step size h; and 
%relative tolerance rTol

%Set up starting values
    t = tspan(1) ;
    tt(:,1) = t ;
    w(:,1) = y0 ;
    ii = 1 ;
    hIn = h ;
% Constants calculated a single time
    c30 = 3/8;
    c31 = 3 / 32 ;
    c32 = 9 / 32 ;
    c40 = 12 / 13 ;
    c41 = 1932 / 2197 ;
    c42 = -7200 / 2197 ;
    c43 = 7296 / 2197 ;
    c51 = 439 / 216 ;
    c52 = -8 ;
    c53 = 3680 / 513 ;
    c54 = -845 / 4104 ;
    c61 = -8 / 27 ;
    c62 = 2 ;
    c63 = -3544 / 2565 ;
    c64 = 1859 / 4104 ;
    c65 = -11 / 40 ;
    cw1 = 25 / 216 ;
    cw3 = 1408 / 2565 ;
    cw4 = 2197 / 4104 ;
    cw5 = -1 / 5 ;
    cz1 = 16 / 135 ;
    cz3 = 6656 / 12825 ;
    cz4 = 28561 / 56430 ;
    cz5 = -9 / 50 ;
    cz6 = 2 / 55 ;
    ce1 = 1 / 360 ;
    ce3 = -128 / 4275 ;
    ce4 = -2197 / 75240 ;
    ce5 = 1 / 50 ;
    ce6 = 2 / 55 ;

    while t+h < tspan(2) %Run until t is to upper bound
        % Calculate RKF values
        s1 = fun(t,w(:,ii));
        s2 = fun(t + 0.25 * h, w(:,ii) + 0.25 * h .* s1);
        s3 = fun(t + c30 * h, w(:,ii) + c31 * h .* s1 + c32 * h .* s2);
        s4 = fun(t + c40 * h, w(:,ii) + c41 * h .* s1 + c42 * h .* s2 + c43 * h .* s3);
        s5 = fun(t + h, w(:,ii) + c51 * h .* s1 + c52 * h .* s2 + c53 * h .* s3 + c54 * h .* s4);
        s6 = fun(t + 0.5 * h, w(:,ii) + c61 * h .* s1 + c62 * h .* s2 + c63 * h .* s3 + c64 * h .* s4 + c65 * h .* s5);
        w(:,ii+1) = w(:,ii) + h .* (cw1 * s1 + cw3 * s3 + cw4 * s4 + cw5 * s5);
        z = w(:,ii) + h .* (cz1 * s1 + cz3 * s3 + cz4 * s4 + cz5 * s5 + cz6 * s6);
        e = h .* norm(ce1 * s1 + ce3 * s3 + ce4 * s4 + ce5 * s5 + ce6 * s6);
        re(ii) = e/norm(w(:,ii+1)) ; %Relative error
        if re(ii) <= rTol %If relative error is less than relative tolerance
            w(:,ii+1) = z ; %Replace w with z
            t = t + h ; %Increase time by a step
            tt(:,ii+1) = t ; %Save time 
            h = hIn ; %Reset step size
            ii = ii+1 ; %Increase indexing
        else %If relative error is too large
            h = h .* .8 * ((rTol*norm(w(:,ii+1)))/e).^.2 ; %Decrease step size and retry
        end
    end
    y = w ;
end

function [ pos ] = lorenz( t,pos )
%lorenz equations. t is a number, pos is a row vector of [ x,y,z]
%values, param is [ sigma , rho , beta ] constants for lorenz
    param = [ 10 , 30 , 8/3 ] ;
    pos = [ param(1) * ( pos(2) - pos(1) ) ; pos(1) * ( param(2) - pos(3) ) - pos(2) ; pos(1) * pos(2) - param(3) * pos(3) ] ;
end