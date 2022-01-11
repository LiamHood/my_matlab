z = [ .1 ; .2 ; .3 ; .4 ; .5 ; .55 ; .6 ; .65 ; .7 ; .75 ; .8 ; .85 ] ; % Vertical location m 
P = [ 100.5 ; 95.6 ; 92.5 ; 93.6 ; 98.4 ; 99.8 ; 98.9 ; 100 ; 99.1 ; 93.0 ; 95.4 ; 99.2 ] ; % Static Pressure kPa 
errP = [ .2 ; .2 ; .4 ; .2 ; 1.3 ; 2.2 ; 1.9 ; 2.4 ; 3.7 ; 2.7 ; .9 ; .3 ] ; % Error kPa
Pstag = 101.3 ; % Stagnation Pressure kPa
rho = 1.225 ; % density in kg/m^3

v = sqrt( ( 2*((Pstag-P)*1000) )/ rho ) ; % Velocity m/s
errv = (errP./P).*v ; % Error for velocity

hold on
errorbar( v , z , errv , 'horizontal' , '.r' ) % Error bar and data plot
plot( v , z , '-b' ) % show lines to aid in seeing trend
title( 'Velocity Profile' )
xlabel( 'Velocity (m/s)' )
ylabel( 'Height from floor (m)' )
hold off

disp( 'I used Bernoulli''s Equation to find velocity which uses the assumptions ' )
disp( 'that the flow is steady, there are no viscous effects, no pumps or turbines' )
disp( 'flow is incompressible, there is no heat transfer, and there are no vortices' )
disp( 'The assumption that there are no vortices and that the flow is steady seem ' )
disp( 'dubious because there seems to be an object in the test section which would ' )
disp( 'likely cause vortices and instability.' )
disp( ' ' )
disp( 'This test could be showing an object in the test section near the top. This' )
disp( 'would explain why the air at the bottom section is fast for longer than the' )
disp( 'air at the top. The data would be clearer if there were more pressure')
disp( 'measurements in the area where the air slows down. As it is it is hard to' )
disp( 'tell what is happening around the increase in speed at .6 m.' )
