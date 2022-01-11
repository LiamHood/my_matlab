clear ; close all ; clc ;
data = readmatrix( 'IMUout.txt' ) ;
time = data( : , 1 )*1e-6 ;
accf = ( data( : , 2 ) - .05 )*98.1 ;
accu = data( : , 3 )*98.1 ;
vela = data( : , 4 ) ;
durationf = data( : , 5 ) ;
durationu = data( : , 6 ) ;
distance = data( : , 7 ) ;
veld = data( : , 8 ) ;

vela = 0 ;
for ii = 2:length( accf )
    vela(ii) = ( accf(ii,1)*( time(ii) - time(ii-1) ) ) + vela(ii-1) ;
end

posa = distance(1) ;
for ii = 2:length( vela )
    posa(ii) = vela(ii)*( time(ii) - time(ii-1) ) + posa(ii-1) ;
end

veld = 0 ;
for ii = 2:length( distance )
    veld(ii) = -( distance(ii) - distance(ii-1) )/( ( time(ii) - time(ii-1) ) ) ;
end

accd = 0 ;
for ii = 2:length( distance )
    accd(ii) = -( veld(ii) - veld(ii-1) )/( ( time(ii) - time(ii-1) ) ) ;
end

figure
hold on
plot( time , accu )
plot( time , accf )
plot( time , accd )
hold off
xlabel( 'Time [s]' )
ylabel( 'Acceleration [cm/s^2]' )
legend( 'Unfiltered Accelerometer' , 'Filtered Accelerometer' , 'Ultrasonic Sensor' )

figure
subplot(2,1,1)
plot( time , vela ) 
xlabel( 'Time [s]' )
ylabel( 'Velocity [cm/s]' )
title( 'Accelerometer' )
subplot(2,1,2)
plot( time , veld )
xlabel( 'Time [s]' )
ylabel( 'Velocity [cm/s]' )
title( 'Ultrasonic Sensor' )

figure
plot( time , distance , time , posa )
xlabel( 'Time [s]' )
ylabel( 'Distance [cm]' )
legend( 'UltraSonic Sensor' , 'Accelerometer' )

