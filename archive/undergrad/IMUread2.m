clear ; close all ; clc ;
data = readmatrix( 'IMUout.txt' ) ;
time = data( 1:426 , 1 )*1e-3 ;
acc = data( 1:426 , 2 ) ;
accUn = data( 1:426 , 3 ) ;
vel = data( 1:426 , 4 ) ;
duration = data( 1:426 , 5 ) ;
durationUn = data( 1:426 , 6 ) ;
distance = data( 1:426 , 7 ) ;
vfd = data( 1:426 , 8 ) ;

dt = time(1) ;
for ii = 2:length(time)
    dt(ii) = time(ii)-time(ii-1) ;
end

vfa = acc(1)*dt(1) ;
for ii = 2:length(acc)
    vfa(ii) = vfa(ii-1) + acc(ii)*dt(ii) ;
end

vfau = accUn(1)*dt(1) ;
for ii = 2:length(accUn)
    vfau(ii) = vfau(ii-1) + accUn(ii)*dt(ii) ;
end

for ii = 10:( length( duration ) )
    duration( ii , 1 ) = mean( duration( (ii-9):(ii) , 1 ) ) ;
end

figure
plot( time , acc , time , accUn )

figure
plot( time , vfd , time , vfa , time , vfau ) 
legend( 'From distance' , 'ReFrom Acceleration' , 'From Acceleration Unfiltered' )

figure
hold on
plot( time , duration )
plot( time , durationUn )
hold off

figure
plot( time , distance ) 
