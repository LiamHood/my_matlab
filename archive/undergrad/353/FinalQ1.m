%% Final Question 1
% Liam Hood
clear; close all; clc;

% Given
abs = .40 ;
emit = .75 ;
nsa = 3.7 ; % surface area of surface normal to solar flux
tsa = 25 ; % total surface area

% changes
cum_d = [ .06 , .10 , .13 , .15 , .17 ] ; % cumulative change 
rabs = abs + cum_d ; % absorption factor
sbc = 5.67e-8 ; % stefan-boltzman constant
Se = 1366 ; 

qsab = rabs .* Se .* nsa ; % heat absorbed
T = ( qsab ./ ( emit * sbc * tsa ) ).^(1/4) ; % Temp of s/c

figure
plot( 1:length(T) , T , 'ob')
title( 'Temperature over time' )
xlabel( 'Years' )
ylabel( 'Spacecraft Temperature (K)' )

disp( 'The spacecraft increases in temperature as the absorption increases' )
disp( 'because the emittance stays the same. This means more thermal energy' )
disp( 'is entering the s/c while the energy leaving remains the same.')

