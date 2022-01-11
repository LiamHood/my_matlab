clear
clc
close all

axis square, xlim([-1 1]), ylim([-1 1]), zlim([-.5 1])
view(125, 50);
grid on

light('position', [1,1,0], 'HandleVisibility', 'off');
lighting gouraud


%% Ready to start?
uiwait(msgbox('Ready to get going?  Remember, the GREEN vector is the initial angular velocity as measured in the body frame.  The YELLOW vector is the angular velocity in the inertial frame.  Don''t forget the angular velocity is constent in the body frame (but not the inertial - our point of view.)','Prolate Case - Stick','modal'));

%% Case 1: Prolate Spinner
title('Case 1: Prolate Spinner')
I_t = 70;
I_a = 40;
w_0 = [0.15 .1 0.25]';  %initial angular velocity in body frame (rad/s)
M_0 = [0; 0; 0];
tspan = 0:.05:20;

Mprolate = simGyroMotion(I_t, I_a, w_0, M_0, tspan)

%% Message Box to continue
uiwait(msgbox('Now on to the OBLATE case.  Same colors apply.','Oblate Case - Frisbee','modal'));
cla

%% Case 2: Oblate Spinner
title('Case 2: Oblate Spinner')
I_t = 40;
I_a = 70;
w_0 = [0.15 .1 0.25]';  %initial angular velocity in body frame (rad/s)
M_0 = [0; 0; 0];
tspan = 0:.05:30;

Moblate = simGyroMotion(I_t, I_a, w_0, M_0, tspan)