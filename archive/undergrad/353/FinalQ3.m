%% Final Question 3
% Liam Hood
clear; close all; clc;

% Given
load( 'Ionosphere' )
alt = [ 2 , 4 , 6 , 10 , 20 ]*100 ;
f = linspace( 10 , 1e4 , 1e3 ) ;
c = 2.99729e8 ;

for ii = 1:5
    denin(ii) = find( Ionosphere == alt(ii) ) ;
end
for ii = 1:5
    for jj = 1:denin(ii)
        TECi(jj) = Ionosphere(jj,2)*200 ;
    end
    TEC(ii) = sum(TECi) ;
    TECi = 0 ;
end
for ii = 1:5
    dt(ii,:) = ( 40.31 * TEC(ii) ) ./ ( c.*f.^2 ) ;
    dr(ii,:) = c*dt(ii,:) ;
end

figure

semilogy( f , dr(1,:) , f , dr(2,:) , f , dr(3,:) , f , dr(4,:) , f , dr(5,:))
title( 'Excess Range vs Frequency' )
xlabel( 'Frequency (Hz)' )
ylabel( 'Excess Range (m)' )
legend( 'Altitude = 200' , '400' , '600' , '1000' , '2000' )
hold off

disp( 'If bending was not neglected the path length would be much longer ')
disp( 'The effect would be more pronounced at the lower frequencies ')
