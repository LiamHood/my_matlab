%% HW 5
% Aero 300
% Liam Hood
clear
close all


%% 6.1.1.e
%Given inputs for function
yprime = @( t , y ) 1/(y^2) ; %function
inter = [ 0 , 1 ] ; %Time interval [ lower bound , upper bound ]
h = .1 ; %Time step
y0 = 1 ; %Intitial y value

[ step , t , w ] = eulerIVP( yprime , inter , h , y0 ) ; %Use Euler method
yAct = @(t) sqrt( 3*t + 1 ) ; %Analytical solution to original function
err = abs(w - yAct(t)) ; %Error of estimate from analytical solution
T_euler = table( step' , t' , w' , yAct(t)' , err' ) ; %Turning data into table
T_euler.Properties.VariableNames = { 'Step' 't_i' 'w_i' 'y_i' 'e_i' } ; %Column heading for table
disp( '6.1.1.e' ) %
disp( T_euler )

%% 6.2.1
%a
    %Given inputs 
    yprime = @(t,y) t ; %function
    inter = [ 0 1 ] ; %time interval
    h = .1 ; %step size
    y0 = 1 ; %intitial y value

    [ step , t , w ] = trapIVP( yprime , inter , h , y0 ) ; %Use ETM
    yAct = @(t) .5*t.^2 + 1 ; %analytic solution to function
    err = abs(w - yAct(t)) ; %error of estimate from true solution
    T_trapa = table( step' , t' , w' , yAct(t)' , err' ) ; %loading data into table
    T_trapa.Properties.VariableNames = { 'Step' 't_i' 'w_i' 'y_i' 'e_i' } ; %table column headings
    disp( '6.2.1.a' )
    disp( T_trapa )
    
%b
    %given inputs 
    yprime = @(t,y) t^2*y ; %given function

    [ step , t , w ] = trapIVP( yprime , inter , h , y0 ) ; %Use ETM
    yAct = @(t)  exp( ( t.^3 ) / 3 ); %analytic solution
    err = abs(w - yAct(t)) ; %error of estimate from true solution
    T_trapb = table( step' , t' , w' , yAct(t)' , err' ) ; %loading data into table
    T_trapb.Properties.VariableNames = { 'Step' 't_i' 'w_i' 'y_i' 'e_i' } ; %table column headings 
    disp( '6.2.1.b' )
    disp( T_trapb )
    
%c
    %given inputs 
    yprime = @(t,y) 2*( t+1 )*y ; %function

    [ step , t , w ] = trapIVP( yprime , inter , h , y0 ) ; %Use ETM
    yAct = @(t) ( 1/exp(1) ) * exp( (t+1).^2 ) ; %analytic solution
    err = abs(w - yAct(t)) ; %error of estimate from true solution
    T_trapc = table( step' , t' , w' , yAct(t)' , err' ) ; %loading data into table
    T_trapc.Properties.VariableNames = { 'Step' 't_i' 'w_i' 'y_i' 'e_i' } ; %table column headings 
    disp( '6.2.1.c' )
    disp( T_trapc )
    
%% 6.2.6
%a large
    %Given inputs 
    yprime = @(t,y) 1 - y^2 ; %function for whole exercise
    inter = [ 0 1 ] ; %interval for whole exercise
    h = .1 ; %large step size
    y0 = 0 ; %starting value a
    
    [ step , t_alarge , w_alarge ] = trapIVP( yprime , inter , h , y0 ) ; %Use ETM
    
%a small 
    h = .05 ; %small step size
    
    [ step , t_asmall , w_asmall ] = trapIVP( yprime , inter , h , y0 ) ; %Use ETM
    
%a exact
    y_a = @(t) tanh( t + atanh( 0 ) ) ; %analytic solution to function with starting value a
    
%b large
    h = .1 ; %large step size
    y0 = -.5 ; %starting value b
    
    [ step , t_blarge , w_blarge ] = trapIVP( yprime , inter , h , y0 ) ; %Use ETM
    
%b small
    %Inputs 
    h = .05 ;
    
    [ step , t_bsmall , w_bsmall ] = trapIVP( yprime , inter , h , y0 ) ; %Use ETM
    
%b exact
    y_b = @(t) tanh( t + atanh( -.5 ) ) ; %analytic solution to function using starting value b
    
%Plots of all estimates and analytic solution
hold on
    plot( t_alarge , w_alarge , '.b' , t_asmall , w_asmall , '.r' , inter(1):.01:inter(2) , y_a( inter(1):.01:inter(2) ) , '-k' ) 
    plot( t_blarge , w_blarge , '.b' , t_bsmall , w_bsmall , '.r' , inter(1):.01:inter(2) , y_b( inter(1):.01:inter(2) ) , '-k' )
    %label graph features
    legend( 'h = .1   y0 = 0' , 'h = .05   y0 = 0' , 'exact   y0 = 0' , 'h = .1   y0 = -.5' , 'h = .05   y0 = -.5' , 'exact   y0 = -.5' , 'Location' , 'southeast' )
    title( 'y'' = 1 - y^2 ' )
    xlabel( 't' )
    ylabel( 'y' )
hold off

%% 6.4.3
%Given inputs
y0 = 1 ; %starting value
inter = [ 0 1 ] ; %interval
h = [ .10 , .05 , .025 ] ; %step sizes

%given functions a through f
yprime4{1} = @( t , y ) t ;
yprime4{2} = @( t , y ) t^2 * y ;
yprime4{3} = @( t , y ) 2 * ( t+1 ) * y ;
yprime4{4} = @( t , y ) 5 * t^4 * y ;
yprime4{5} = @( t , y ) y^-2 ;
yprime4{6} = @( t , y )  t^3 / y^2 ;

%analytic solution to a through f
y{1} = @(t) .5 * t.^2 + 1 ;
y{2} = @(t) exp( ( t.^3) / 3 ) ;
y{3} = @(t) exp( ( t+1 ).^2 -1 ) ;
y{4} = @(t) exp( t.^5 ) ;
y{5} = @(t) ( 3*t + 1 ).^(1/3) ;
y{6} = @(t) ( .75*t.^4 + 1 ).^(1/3) ;

for jj = 1:6 %Runs for all given functions a-f
    yprime = yprime4{jj} ; %assign appropriate function handle to be referred to by rk4
    figure
    hold on
    for kk = 1:3 %Evaluates each of the functions using each of the three step sizes
        [ step4{kk} , t4{kk} , w4{kk} ] = rk4( yprime , inter , h(kk) , y0 ); %Uses RK4
        plot( t4{kk} , w4{kk} ) %Plots estimate
        
        % Create a table for each step size and each function. Tables are
        % organized together as a cell array. Each row corresponds to one
        % of the given functions with row 1 being function 1, row 2 being
        % function 2, and so on. Each column corresponds to a step size
        % decreasing from left to right
        err = abs(w4{kk} - y{jj}(t4{kk})) ; %error of estimate from true solution
        T_rk4{jj,kk} = table( step4{kk}' , t4{kk}' , w4{kk}' , y{jj}(t4{kk})' , err' ) ; %loading data into a table
        T_rk4{jj,kk}.Properties.VariableNames = { 'Step' 't_i' 'w_i' 'y_i' 'e_i' } ; %table column headings
        disp( '6.4.3' )
        disp([ 'RK4 approximation of ' , func2str( yprime ) ,'with step size ' , num2str( h(kk) ) ]) 
        disp( T_rk4{jj,kk} )
    end
    t_int = linspace( inter(1) , inter(2) , 100 ) ; %points to evaluate analytic solution on
    plot( t_int , y{jj}( t_int ) ) %plot true/analytic solution
    %label plot features
    title([ 'y''' , func2str( yprime ) ])
    xlabel( 't' )
    ylabel( 'y' )
    legend( 'h = .1' , 'h = .05' , 'h = .025' , 'actual function' )
    hold off
end
    
%% Functions

function [ step , t , w ] = eulerIVP( yprime , inter , h , y0 )
% Uses Euler''s Method to evaluate IVPs. yprime is a function handle for
% the given differiental. inter is the interval from low to high as a
% vector. h is step size. y0 is the initial y value. The outputs are w
% values representing Euler estimates of y values with corresponding step
% number and time value

    w(1) = y0 ; %Start estimate with the given initial value
    t(1) = inter(1) ; %Start time at lower bound of interval
    ii = 1 ;%Start indexing variable at 1
    
    while t+h < inter(2) %Runs as long as t is less than the upper limit of the interval
        step(ii+1) = ii ; %steps begins at 0 because no steps are needed to find the initial value. Increases as index does 
        w(ii+1) = w(ii) + h*yprime( t(ii) , w(ii) ) ; %Euler Method formula
        t(ii+1) = t(ii) + h ; %Time increases by step size each step
        ii = ii + 1 ; %Index increases by one each iteration
    end

end

function [ step , t , w ] = trapIVP( yprime , inter , h , y0 )
% Uses Explicit Trapezoid Method to evaluate IVPs. yprime is a function handle for
% the given differiental. inter is the interval from low to high as a
% vector. h is step size. y0 is the initial y value. The outputs are w
% values representing ETM estimates of y values with corresponding step
% number and time value

    w(1) = y0 ; %Start estimate with the given initial value
    t(1) = inter(1) ; %Start time at lower bound of interval
    ii = 1 ; %Start indexing variable at 1
    
    while t+h < inter(2) %Runs as long as t is less than the upper limit of the interval
        step(ii+1) = ii ; %steps begins at 0 because no steps are needed to find the initial value. Increases as index does 
        w(ii+1) = w(ii) + .5*h*( yprime( t(ii) , w(ii) ) + yprime( t(ii)+h , w(ii) + h*yprime( t(ii) , w(ii) )))  ; %ETM formula
        t(ii+1) = t(ii) + h ; %Time increases by step size each step
        ii = ii + 1 ; %Index increases by one each iteration
    end

end

function [ step , t , w ] = rk4( yprime , inter , h , y0 )
% Runge Katta fourth order IVP solver.  yprime is a function handle for
% the given differiental. inter is the interval from low to high as a
% vector. h is step size. y0 is the initial y value. The outputs are w
% values representing ETM estimates of y values with corresponding step
% number and time value
    w(1) = y0 ; %Start estimate with the given initial value
    t(1) = inter(1) ; %Start time at lower bound of interval
    ii = 1 ; %Start indexing variable at 1
    
    while t < inter(2) %Runs as long as t is less than the upper limit of the interval
        step( ii+1 ) = ii ; %steps begins at 0 because no steps are needed to find the initial value. Increases as index does 
            %slope estimates withing each step. Used to estimate next step
            s1 = yprime( t(ii) , w(ii) ) ;
            s2 = yprime( t(ii) + (h/2) , w(ii) + (h/2)*s1 ) ;
            s3 = yprime( t(ii) + (h/2) , w(ii) + (h/2)*s2 ) ;
            s4 = yprime( t(ii) + h , w(ii) + h*s3 ) ;
        w( ii+1 ) = w(ii) + ( h/6 ) * ( s1 + 2*s2 + 2*s3 + s4 ) ; % 4th order Runge-Katta formula
        t( ii+1 ) = t(ii) + h ; %Time increases by step size each step
        ii = ii + 1 ; %Index increases by one each iteration
    end
end