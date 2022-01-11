%% HW1 Aero 402
% Liam Hood

clear ; close all ; clc ;

type = [ "Cold Gas" , "Chemical" , "Nuclear" , "Electric" ] ;
Isp = [ 100 , 400 , 1000 , 5000 ] ; % [s]
finert = linspace( 0 , 1 , 1e3 ) ;

for ii = 1:4 
    deltav(ii,:) = Isp(ii)*log( 1./finert )*9.81 ;
end

figure
semilogy( finert , deltav(1,:) , 'b' , finert , deltav(2,:) , 'r' , finert , deltav(3,:) , 'm' , finert , deltav(4,:) , 'g' )
hold on

area( 0:1 , 30520*[1,1] , 'FaceColor' , [ .6 .6 0 ] ) % there and back again neptune
area( 0:1 , 12990*[1,1] , 'FaceColor' , [ .6 .4 0 ] ) % leo-mercury
area( 0:1 , [ 3940 , 3940 ]  , 'FaceColor' , [ .6 .2 0 ] ) % leo-moon
area( 0:1 ,[ 50 , 50 ]  , 'FaceColor' , [ .6 0 0 ] ) % attitude
semilogy( finert , deltav(1,:) , 'b' , finert , deltav(2,:) , 'r' , finert , deltav(3,:) , 'm' , finert , deltav(4,:) , 'g' )
xlabel( 'Inert Mass Fraction' )
ylabel( 'Delta V [m/s]' )
legend( type  )
grid
hold off