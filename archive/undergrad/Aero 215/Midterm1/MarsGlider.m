%Liam Hood
%Aero 215 Aircraft Midterm
%Lift and Drag of Martian glider
    %Calculates lift, drag, and L/D for a glider with given coefficients of
    %lift and drag, as well as wing area and starting altitude and velocity
    clear all;
    clc;   
    %%Birdy Glider
        disp('Birdy')
        %Input glider characteristics
        CL = .65; %Coefficient of Lift
        CD = .02; %Coefficient of Drag
        S = 48;   %m^2
        v = 45;   %velocity at launch in m/s
        h = 16000;%height of launch in m
        MarsGliderFunction( CL , CD , S , v , h );%runs function to calculate lift, drag, and L/D

    %%Raptor Glider
        disp('Raptor')
        %Input glider characteristics
        CL = .80; %Coefficient of Lift
        CD = .012; %Coefficient of Drag
        S = 110;   %m^2
        v = 35;   %velocity at launch in m/s
        h = 16000;%height of launch in m
        MarsGliderFunction( CL , CD , S , v , h );%runs function to calculate lift, drag, and L/D

    %%Boomer Glider
        disp('Boomer')
        %Input glider characteristics
        CL = 1.35; %Coefficient of Lift
        CD = .023; %Coefficient of Drag
        S = 202;   %m^2
        v = 18;   %velocity at launch in m/s
        h = 16000;%height of launch in m
        MarsGliderFunction( CL , CD , S , v , h );%runs function to calculate lift, drag, and L/D    





