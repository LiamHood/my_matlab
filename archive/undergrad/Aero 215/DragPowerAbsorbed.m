%% Aero 215 HW2
% Liam Hood
% 10/13/17
% Drag & Power Absorbed
close all; clear all; clc
%format short
%% Initial Setup
v_max = 300*1.69 ; % [ft/s] Maximum Incoming Velocity (300 knots)
inputs.WS = 70 ; % [lbf/ft^2] Wing Loading
inputs.AR = 8 ; % Aspect Ratio
inputs.Cdp = 0.02 ; % Parasite Drag Coefficient, Clean
inputs.e = 0.8 ; % Oswald Efficiency
inputs.cl_max = 1.2 ; % Maximum Lift Coefficient, Clean
inputs.S = 200 ; % [ft^2] Wing Reference Area
inputs.W = inputs.WS*200 ; % [lbf] Weight
inputs.rho = 0.00238 ; % [slug/ft^3] Density @ Sea Lvl
%% Sea Level, Clean
disp( 'Sea Level, Clean' )
[v, v_sl, D_sl, P_sl] = HW2_Drag_Power_Hood_Liam(v_max, inputs);


%% Cruise Altitude, Clean (10,000 ft)
disp( 'Cruise Altitude, Clean (10,000 ft)' )
%reset inputs
inputs.rho = .001756; % [slug/ft^3]
[vca, v_ca, D_ca, P_ca] = HW2_Drag_Power_Hood_Liam(v_max, inputs);


%% Sea Level, Dirty
disp( 'Sea Level, Dirty' )
%reset inputs
inputs.rho = .00238; % [slug/ft^3]
inputs.Cdp = .03; % Parasite Drag Coefficient
insputs.cl_max = 2; % Maximum Lift Coefficient
[vdl, v_dl, D_dl, P_dl] = HW2_Drag_Power_Hood_Liam(v_max, inputs);


%% Sea Level, Varied Wing Loading
disp( 'Sea Level, Varied Wing Loading' )
%reset inputs
inputs.Cdp = .02; % Parasite Drag Coefficient, Clean
inputs.cl_max = 1.2; % Maximum Lift Coefficient, Clean
inputs.WS = 75; % [lbf/ft^2] Wing Loading
inputs.W = inputs.WS*200; % [lbf] Weight
[vvwl, v_vwl, D_vwl, P_vwl] = HW2_Drag_Power_Hood_Liam(v_max, inputs);


%% Plotting
subplot( 1 , 2 , 1) %First of the two subplots is drag in all flight conditions
plot( v , D_sl , 'g' , vca , D_ca , vdl , D_dl , vvwl , D_vwl , ':b' )
%labels graph
legend ( 'Sea Level' , 'Cruise Altitude' , 'Landing Configuration' , 'Higher wing loading' )
xlabel( 'Velocity in ft/s' )
ylabel( 'Drag Force in lbf' )
title( 'Drag vs. Velocity' )

subplot( 1 , 2 , 2 ) %Second subplot is power in all flight conditions
plot( v , P_sl , 'g' , vca , P_ca , vdl , P_dl , vvwl , P_vwl , ':b' )
%labels graph
legend ( 'Sea Level' , 'Cruise Altitude' , 'Landing Configuration' , 'Higher wing loading' )
xlabel( 'Velocity in ft/s' )
ylabel( 'Power Required in kW' )
title( 'Power vs. Velocity' )

