%Liam Hood Aero 215 Section 3 Homework 1
%Stand Atmosphere (1959) Model
close all;
clear all;
clc;

h = 0; %h is sltitude in meters less than 100000
[ T , P , rho ] = stdatm_HOOD_LIAM( h ); %Calculates temperature, pressure, and density at any altitude below 100km
Td = [ 'The temperature at '  num2str(h)  ' meters  is  '  num2str(T(1))  ' K' ] ; %Adds context to temperature
disp ( Td ) 
Pd = [ 'The pressure at ' , num2str(h) , ' meters  is  ' , num2str(P(1)) , ' kPa' ]; %Adds context to pressure
disp( Pd )
rhod = [ 'The density at  ' , num2str(h) , ' meters  is  ' , num2str(rho(1)) , ' kg/m^3' ]; %Adds context to density
disp ( rhod )