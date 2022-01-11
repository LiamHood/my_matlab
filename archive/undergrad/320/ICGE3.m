
clear ; close all ; clc ;
L = 1 ;
p = .25 ;
m = 1 ;
B = .5 ;
g = 9.81 ;
x0 = [ 0 , .5 ] ;
opts = [ 'AbsTol' , 1e-8 , 'RelTol' , 1e-8 ] ;
tspan = [ 0 , 50 ] ;
[ t6 , x6 ] = ode45( @DoublePendulum , tspan , x0 , opts , L , p , m , B , g ) ;
    % Plots 
        figure
        plot( t6 , x6(:,1) )
        title( 'Theta vs Time' )
        xlabel( 'Time (s)' )
        ylabel( 'Theta (radians)' )

        figure
        plot( t6 , x6(:,2) )
        title( 'Theta Dot vs Time' )
        xlabel( 'Time (s)' )
        ylabel( 'Theta Dot (radians/s)' )

        figure 
        plot( x6(:,1) , x6(:,2) )
        title( 'Theta Dot vs Theta ' )
        xlabel( 'Theta (radians)' )
        ylabel( 'Theta Dot (radians/s)' )
        
 p = -.25:.125:.25 ;    
       
    [ t1 , x1 ] = ode45( @DoublePendulum , tspan , x0 , opts , L , p(1) , m , B , g ) ;
    [ t2 , x2 ] = ode45( @DoublePendulum , tspan , x0 , opts , L , p(2) , m , B , g ) ;
    [ t3 , x3 ] = ode45( @DoublePendulum , tspan , x0 , opts , L , p(3) , m , B , g ) ;
    [ t4 , x4 ] = ode45( @DoublePendulum , tspan , x0 , opts , L , p(4) , m , B , g ) ;
    [ t5 , x5 ] = ode45( @DoublePendulum , tspan , x0 , opts , L , p(5) , m , B , g ) ;
    % Plots 
        figure
        hold on
        plot( t1 , x1(:,1) )
        plot( t2 , x2(:,1) )
        plot( t3 , x3(:,1) )
        plot( t4 , x4(:,1) )
        plot( t5 , x5(:,1) )
        title( 'Theta vs Time' )
        xlabel( 'Time (s)' )
        ylabel( 'Theta (radians)' )
        legend( 'p = -.25' , 'p = -.125' , 'p = 0' , 'p = .125' , 'p = .25' )
        hold off

        figure
        hold on
        plot( t1 , x1(:,2) )
        plot( t2 , x2(:,2) )
        plot( t3 , x3(:,2) )
        plot( t4 , x4(:,2) )
        plot( t5 , x5(:,2) )
        title( 'Theta Dot vs Time' )
        xlabel( 'Time (s)' )
        ylabel( 'Theta Dot (radians/s)' )
        legend( 'p = -.25' , 'p = -.125' , 'p = 0' , 'p = .125' , 'p = .25' )
        hold off

        figure 
        hold on 
        plot( x1(:,1) , x1(:,2) )
        plot( x2(:,1) , x2(:,2) )
        plot( x3(:,1) , x3(:,2) )
        plot( x4(:,1) , x4(:,2) )
        plot( x5(:,1) , x5(:,2) )
        title( 'Theta Dot vs Theta' )
        xlabel( 'Theta (radians)' )
        ylabel( 'Theta Dot (radians/s)' )
        legend( 'p = -.25' , 'p = -.125' , 'p = 0' , 'p = .125' , 'p = .25' )
        hold off

%% Functions

function [ dx ] = DoublePendulum( t , x , L , p , m , B , g )
    dx = [ x(2) ; -(B/(2*m*L^2))*x(2) - (g*p/L^2)*sin(x(1)) ] ;
end

