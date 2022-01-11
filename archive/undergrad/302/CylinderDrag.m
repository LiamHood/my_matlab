clear ; clc ; close all ;
load( 'Cy_500_1.mat' )
P_scanivalve = P ;
n = 24 ;
% turn pressure into logically numbered matrix
for ii = 1:n
    P(:,ii) = P_scanivalve(:,ii+2) ; 
end

% Set up 
d2r = 180/pi ; % degrees to radians
sep = 15*d2r ; % seperation angle in radians
r = .15/2 ; % radius in meter
% l = 1.22 ; % length
theta = sep ; % Initial 

% Find pressure force at each port
for jj = 1:n
    theta(jj) = theta(1) + sep*(jj-1) ; % angle of port
    fx(:,jj) = - P(:,jj) * r * sep * cos( theta(jj) ) ;
end

% Find mean, standard deviation, and standard error
for kk = 1:n
    fx_mean(kk) = mean( fx(:,kk) ) ;
    dev(kk) = sum( fx(:,kk)-fx_mean(kk) ).^2 ;
    sig(kk) = sqrt( dev(kk) / ( n - 1 ) ) ;
    alpha(kk) = sig(kk)/sqrt(n) ;
end

figure
hold on
title( 'Drag force at each pressure port' )
xlabel( 'Angle of pressure port clockwise from trailing edge horizontal (degrees)' )
ylabel( 'Drag force per unit width (N/m)' )
errorbar( theta/d2r , fx_mean , 3*alpha , '.-b' ) 
plot( theta/d2r , fx_mean , '.r')
hold off

drag = sum( fx_mean ) ;
err = sqrt( sum( alpha.^2 ) ) ;

disp([ 'Drag is ' , num2str( drag ) , ' N/m with a standard error of ' , num2str( err ) ])

