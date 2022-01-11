clear ; close all ;
tspan = [ 0 , 25 ] ;
x0 = [ 1 ; 0 ] ;
nf = 1 ;
dr = [ 0 , .25 , 1 , 2 ] ;
part = [ "undamped" , "underdamped" , "critically damped" , "overdamped" ] ;


for ii = 1:4 
    options = odeset( 'AbsTol' , 1e-8 , 'RelTol' , 1e-8 ) ;
    [ t , x ] = ode45( @SpringDamper , tspan , x0 , options , nf , dr(ii) ) ;
    figure 
    plot( t , x(:,1) )
    title([ 'Position vs Time ' , part(ii) ])
    xlabel( 'Time' )
    ylabel( 'Position' )
    figure
    plot( t , x(:,2) )
    title([ 'Velocity vs Time ' , part(ii) ])
    xlabel( 'Time' )
    ylabel( 'Velocity' )
    figure
    plot( x(:,1) , x(:,2) )
    title([ 'State-Space ' , part(ii) ])
    xlabel( 'Position' )
    ylabel( 'Velocity' )
end

disp( '11' )
disp( 'The best scenario is is critically damped because it achieves steady ')
disp( 'state soonest but overdamped is better thannunderdamped' )