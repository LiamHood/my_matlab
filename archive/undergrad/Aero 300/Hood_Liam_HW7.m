%% HW 7
% Aero 300
% Liam Hood

clear
close all

%% 8.1.1b
disp( '8.1.1b' )
%Inputs
h = .1 ;
k = .002 ;
xspan = [ 0 1 ] ;
tspan = [ 0 1 ] ;
ubeg = @(x) exp(x) ;
uleft = @(t) exp( 2*t ) ;
uright = @(t) exp( 2*t + 1 ) ;
D = 2 ;

[ w , x , t ]= parabolicFDM( ubeg , uleft , uright , tspan , xspan , h , k , D); %Using FDM
%Graphing solution
mesh( x , t , w )
title( 'Heat Function by FDM' )
xlabel( 'x' )
ylabel( 'time' )
zlabel( 'Temperature' )

% Checking against true solution
uAct = @(x,t) exp( 2*t + x ) ;
s = 1000 ;
tAct = linspace( tspan(1) , tspan(2) , s ) ;
xAct = linspace( xspan(1) , xspan(2) ,  s) ;
for ii = 1:s
    for jj = 1:s
        uAct1(ii,jj) = uAct( xAct(ii) , tAct(jj) ) ;
    end
end
%Graphing true solution
figure
mesh( tAct , xAct , uAct1' ) ;
title( 'Actual Heat Function' )
xlabel( 'x' )
ylabel( 'time' )
zlabel( 'Temperature' )

[ w , x , t ]= parabolicFDM( ubeg , uleft , uright , tspan , xspan , h , .0027 , D); %Using FDM
%Graphing solution
figure
mesh( x , t , w )
title( 'Heat Function by FDM with large time step' )
xlabel( 'x' )
ylabel( 'time' )
zlabel( 'Temperature' )
disp( 'The graph behaves strangely when k is greater than .0026. The solution ' )
disp( 'gets wavy on one side and flat on the other instead of the smooth curve' )
disp( ' it is meant to have.')

%% 8.1.3b
disp( '8.1.3b' )
%Inputs
h = .02 ;
k1 = .02 ;
k2 = .01 ;
k3 = .005 ;
k4 = .002 ; 
xspan = [ 0 1 ] ;
tspan = [ 0 1 ] ;
ubeg = @(x) exp(x) ;
uleft = @(t) exp( 2*t ) ;
uright = @(t) exp( 2*t + 1 ) ;
D = 2 ;

[ w1 , x1 , t1 ]= parabolicBDM( ubeg , uleft , uright , tspan , xspan , h , k1 , D); %Using BDM
[ w2 , x2 , t2 ]= parabolicBDM( ubeg , uleft , uright , tspan , xspan , h , k2 , D); %Using BDM
[ w3 , x3 , t3 ]= parabolicBDM( ubeg , uleft , uright , tspan , xspan , h , k3 , D); %Using BDM

uAct2 = [ uAct(.5,1) ; uAct(.5,1) ; uAct(.5,1) ] ;
k = [ k1 ; k2 ; k3 ] ;
wd = [ w1( length(w1) , 26 ) ; w2( length(w2) , 26 ) ; w3( length(w3) , 26 ) ] ;
e = abs( uAct2 - wd );
data = table( k , uAct2 , wd , e ) ;
data.Properties.VariableNames = { 'Time_step'  'Actual_value'  'Estimated_value'  'Error' } ;
disp( data )

%% 8.2.1b
%Inputs
c = 2 ;
h = .05 ;
k = h/c ;
ubeg = @(x) exp( -x ) ;
ubegt = @(x) -2*exp( -x ) ;
uleft = @(t) exp( -2*t ) ;
uright = @(t) exp( -1 - 2*t ) ;
%Use function
[ ww , xw , tw] = ellipticFDM( ubeg , ubegt ,uleft , uright , tspan , xspan , h , k , D) ; 
%Graph results
figure
mesh( xw , tw , ww )
title( 'Wave function' )
xlabel( 'x' )
ylabel( 'time' )
zlabel( 'amplitude' )

%% Functions

function [ w , x , t ] = parabolicFDM( ubeg , uleft , uright , tspan , xspan , h , k , D)
% Uses Forward difference method on parabolic PDEs
    n = ( tspan(2) - tspan(1) ) / k ; 
    m = ( xspan(2) - xspan(1) ) / h - 1;

    sigma = D*k/(h^2) ; %define sigma
    a = diag( 1 - 2*sigma*ones( m , 1 )) + diag( sigma*ones(m-1,1) , 1 ) + diag( sigma*ones(m-1,1) , -1 ) ; %define a
    left = uleft( tspan(1) + (0:n)*k ) ; %left side boundary
    right = uright( tspan(1) + (0:n)*k) ; %right side boundary
    w(:,1) = ubeg( xspan(1) + (1:m)*h )' ; %initial condition

    for ii = 1:n 
        w(:,ii+1) = a*w(:,ii)+sigma*[left(ii);zeros(m-2,1);right(ii)] ; %Apply formula
    end
    w = [ left ; w ; right ]' ;
    x = ( 0:m+1 )*h ;
    t = ( 0:n )*k ;
end

function [ w , x , t ] = parabolicBDM( ubeg , uleft , uright , tspan , xspan , h , k , D)
% Uses backward difference method on parabolic PDEs
    n = ( tspan(2) - tspan(1) ) / k ; 
    m = ( xspan(2) - xspan(1) ) / h - 1;

    sigma = D*k/(h^2) ; %define sigma
    a = diag( 1 + 2*sigma*ones( m , 1 )) + diag( -sigma*ones(m-1,1) , 1 ) + diag( -sigma*ones(m-1,1) , -1 ) ; %define a
    left = uleft( tspan(1) + (0:n)*k ) ; %left side boundary
    right = uright( tspan(1) + (0:n)*k) ; %right side boundary
    w(:,1) = ubeg( xspan(1) + (1:m)*h )' ; %initial condition

    for ii = 1:n 
        w(:,ii+1) = a*w(:,ii)+sigma*[left(ii);zeros(m-2,1);right(ii)] ; %Apply formula
    end
    w = [ left ; w ; right ]' ;
    x = ( 0:m+1 )*h ;
    t = ( 0:n )*k ;
end

function [ w , x , t ] = ellipticFDM( ubeg , ubegt ,uleft , uright , tspan , xspan , h , k , c )
% Uses finite difference method on parabolic PDEs
    n = ( tspan(2) - tspan(1) ) / k ; 
    m = ( xspan(2) - xspan(1) ) / h - 1;

    sigma = c*k/h ; %define sigma
    a = diag( 2 - 2*sigma^2*ones( m , 1 )) + diag( sigma^2*ones(m-1,1) , 1 ) + diag( sigma^2*ones(m-1,1) , -1 ) ; %define a
    left = uleft( tspan(1) + (0:n)*k ) ; %left side boundary
    right = uright( tspan(1) + (0:n)*k) ; %right side boundary
    w(:,1) = ubeg( xspan(1) + (1:m)*h )' ; %initial condition
    w(:,2) = .5*a*w(:,1) + k*ubegt( xspan(1) + (1:m)*h )' + .5*sigma^2*[left(1);zeros(m-2,1);right(1)] ; %second condition
    
    for ii = 2:n 
        w(:,ii+1) = a*w(:,ii)+sigma*[left(ii);zeros(m-2,1);right(ii)] ; %Apply formula
    end
    w = [ left ; w ; right ]' ;
    x = ( 0:m+1 )*h ;
    t = ( 0:n )*k ;
end