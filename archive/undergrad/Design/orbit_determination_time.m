clear ; close all ; clc ;

minutes = [ 10, 25, 50, 100] ;
R = [.3337, .3430, .3508, .3664] ;
I = [.1235, .1385, .1594, .1942] *5 ;
C = [.0792, .0900, .1093, .1557] ;

figure
title('Two-Sigma Position Uncertainty')
hold on 
plot( minutes , R*2 )
plot( minutes , I*2 ) 
plot( minutes , C*2 )
legend('Radial', 'In-Track', 'Cross-Track', 'Location', 'west')
xlabel('Time Between Observations [minutes]')
ylabel('Position Uncertainty [meters]')
grid on