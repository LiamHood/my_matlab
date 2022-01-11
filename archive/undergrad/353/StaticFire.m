clear; close all; clc;
load( 'rocket.txt' ) ;
time = rocket(:,1) ;
thrust_lb = rocket(:,2) ;
thrust = thrust_lb / .22481 ;
time = time( (length(time)/4):(length(time)/2) ) ;
thrust = thrust( (length(thrust)/4):(length(thrust)/2) ) ;

burn_thrust_i = find( thrust > 1.68 ) ;
burn_thrust = thrust( burn_thrust_i ) ;
burn_time = time( burn_thrust_i ) - time( min( burn_thrust_i ) ) ;

tstep = .002 ; % seconds
dimpulse = burn_thrust * tstep ;
impulse = sum( dimpulse ) ;
f_avg = mean( burn_thrust ) ;



% Adjust down to give similar average thrust and impulse to published data
burn_thrustn = thrust( burn_thrust_i ) - .8 ;

dimpulse = burn_thrustn * tstep ;
impulsen = sum( dimpulse ) ;
f_avgn = mean( burn_thrustn ) ;


estesB6 = [ 0.023 0.688 ; ...
0.057 2.457 ; ...
0.089 4.816 ; ...
0.116 7.274 ; ...
0.148 9.929 ; ...
0.171 12.140 ; ... 
0.191 11.695 ; ...
0.200 10.719 ; ...
0.209 9.240 ; ...
0.230 7.667 ; ...
0.255 6.488 ; ...
0.305 5.505 ; ...
0.375 4.816 ; ...
0.477 4.620 ; ...
0.580 4.620 ; ...
0.671 4.521 ; ...
0.746 4.226 ; ...
0.786 4.325 ; ...
0.802 3.145 ; ...
0.825 1.572 ; ...
0.860 0.00 ] ;

figure
plot( estesB6(:,1) , estesB6(:,2) , burn_time , burn_thrust , burn_time , burn_thrustn )
title( 'Thrust Curves' )
xlabel( 'Thrust (N)' )
ylabel( 'Time (s)' )
legend( 'Published Data' , 'Raw Burn Data' , 'Shifted down to have similar average force' )
disp( 'B6-6' )
disp( 'B means that the impulse is 5 N-sec' )
disp( 'First 6 says which B engine it is' )
disp( 'Second 6 says there is a 6 second delay before the ejection charge' )
thrusttab = [ "Average thrust " , num2str(f_avg) , " Newtons" ] ;
impulsetab = [ "Impulse " , num2str(impulse) , " Newton*seconds" ] ;

massloss = 5.6 ; % grams lost to black powder according to published data. Ours lost 9.7 but this wasn't all due to the thrust
massflow = massloss/max(burn_time) ;
massflowtab = [ "Mass Flow Rate " , num2str(massflow) , " gram*seconds" ] ;

massflowk = massflow/1000 ;
ve = f_avg/massflowk ;
Isp = ve/9.81 ;
vetab = [ "Effective exhaust velocity " , num2str( ve ) , " m/s" ] ;
Isptab = [ "Specific Impulse " , num2str(Isp) , " 1/s" ] ;
table = [ thrusttab ; impulsetab ; massflowtab ; vetab ; Isptab ] ;
disp( table )