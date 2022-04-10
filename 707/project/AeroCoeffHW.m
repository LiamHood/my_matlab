clc;clear all;close all;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')

%% Model geometry (Bixler 2)
W = 2; % lb Weight of the airplane
rho = 0.002378;% density of air
U1 = 45; % ft/s forward velocity
S = 2.6; % ft^s wing area
S_h = 0.49; % 144 in^2=1 ft^2
IYY = 0.02274; % Slug*ft^2
cbar = 0.74; % chord of the wing
x_cg = 0.202; % distance from LE of wing to the c.g location of the aircraft
x_ac = 0.15; % Assumption
x_ac_h = 1.708; % the distance from LE of the wing to AC of Horizon tail
xbar_ac_wf = x_ac/cbar;
xbar_ac_h = x_ac_h/cbar;
xbar_cg = x_cg/cbar;
qbar=1/2*rho*U1^2; % dynamic pressure
g = 32.2;
m = W/g; % slug mass of airplane
l_h = x_ac_h- x_cg; % 3.154 (pg 145)
Vbar_h = S_h/S*(xbar_ac_h - xbar_cg);

%% Assumptions
CMu = 0; % Page 138 figure 3.55
CM1 = 0;
CM0 = -0.074;
eta_h = 0.9;
deda = 0.25;
CL0_wf = 0.35;
CL0_h = 0;
CLa_h = 0.08*180/pi;
CLa_wf = 0.09*180/pi;
epsilon0 = 3*pi/180;
tau_e = 0.4111;
ih = 0;
de = -7.152; % trim elevator
at = 0; % thrust angle
gammae = 0; % flight path
alpha = 0.5; % trim alpha
a = 1115.64; % speed of sound in ft/sec
M1 = U1/a;

%% Lift
CL0 = CL0_wf - CLa_h*eta_h*S_h/S*epsilon0 + CL0_h; % 3.23 (pg 79)
CLa = CLa_wf + CLa_h*eta_h*S_h/S*(1-deda); % 3.24 (pg 79)
CLadot = 2*CLa_h*eta_h*Vbar_h*0.1;% pg 141
CLq = 2*CLa_h*l_h/cbar*eta_h*S_h/S; % 3.153 (pg 145)
CLde = CLa_h*eta_h*S_h/S*tau_e; % 3.26 (pg 79)
CLih = CLa_h*eta_h*S_h/S; % 3.25 (pg 79)
CL = CL0 + CLa*alpha + CLih*ih + CLde*de; % 3.12 (pg 77)
CLu = M1^2/(1-M1^2)*CL;
CL1 = W/qbar/S;

%% Drag
CDu=0; % pg 133 figure 3.53 on pg135

%% Pitching Moment
% CM0 = Cmac_wf + CL0_wf*(xbar_cg - xbar_ac_wf) + CLa_h*eta_h*Vbar_h*epsilon0; % 3.34 (pg 84)
% CMq = -2.2*CLa_h*eta_h*Vbar_h*(xbar_ac_h - xbar_cg); % 3.158 (pg 145)
% CMih = -CLa_h*eta_h*Vbar_h; % 3.36 (pg 85)

CMa = CLa_wf*(xbar_cg - xbar_ac_wf) - CLa_h*eta_h*Vbar_h*(1-deda); % 3.35 (pg 85)
CMde = -CLa_h*eta_h*Vbar_h*tau_e; % 3.37 (pg 85)
CMq = -2*CLa_h*eta_h*Vbar_h*(xbar_ac_h - xbar_cg); % 3.157 (pg 145)
CMadot = -2*CLa_h*eta_h*Vbar_h*(xbar_ac_h - xbar_cg)*0.1; %


%% Rolling Moment
% Clb_wf = 0; % pg 96
% Clb_h = 0; % 3.61 (pg 100)
% Clb_v = -CLa_v*(1-dsdb)*eta_v*S_v*z_v/(S*b); % 3.65 (pg101)
% Clb = Clb_wf + Clb_h + Clb_v; % 3.52 (pg96)
% Clda = 0; % pg 103
% Cldr = CLa_v*alpha_dr*qbar_v*S_v*x_v/(S*b); % 4.71 (pg 108)
% Clp_wf = 0; % pg 152
% Clp_h = 0; % pg 152
% Clp_v = -2*CLa_v*(z_v/b)^2*eta_v*(S_v/S); % 3.180 (pg 152)
% Clp = Clp_wf + Clp_h + Clp_v; % 3.177 (pg 150)
% Clr_wf = 0; % pg 158
% Clr_h = 0; % pg 158
% Clr_v = CLa_v*(2*x_v*z_v/b)*eta_v*S_v/S; % 3.191 (pg 158)
% Clr = Clr_wf + Clr_h + Clr_v; % 3.190 (pg 158)
% 
%% Side Force
% CY0 = 0; % for symmetric (pg 110)
% CYb_w = 0; % pg 111
% CYb_f = 0; % pg 111
% CYb_v = CLa_v*(1-dsdb)*eta_v*S_v/S; % 3.77 (pg 111)
% CYb = CYb_w + CYb_f + CYb_v; % 3.75 (pg 110)
% CYda = 0; % pg 111
% CYdr = CLa_v*alpha_dr*eta_v*S_v/S; % 3.79 (pg 113)
% CYp = -2*CLa_v*(z_v/b)*eta_v*S_v/S; % 3.176 (pg 150)
% CYr_wf = 0; % pg 157
% CYr_v = CLa_v*(2*x_v/b)*eta_v*S_v/S; % 3.189 (pg 157)
% CYr = CYr_wf + CYr_v; % 3.186 (pg 157)
% 
%% Yawing Moment
% Cnb_w = 0; % pg 115
% Cnb_f = 0; % pg 115
% Cnb_v = CLa_v*(1-dsdb)*eta_v*S_v*x_v/(S*b); % 3.85 (pg 111)
% Cnb = Cnb_w + Cnb_f + Cnb_v; % 3.83 (pg 115)
% Cnda = 0; % pg 117
% Cndr = -CLa_v*alpha_dr*eta_v*S_v*x_v/(S*b); % 3.87 (pg 120)
% Cnp_wf = 0; % pg 152
% Cnp_v = 2*CLa_v*(z_v/b)*(x_v/b)*eta_v*S_v/S; % 3.182 (pg 154)
% Cnp = Cnp_wf + Cnp_v; % 3.181 (pg 152)
% Cnr_wf = 0; % pg 160
% Cnr_v = -CLa_v*(2*x_v^2/b)*eta_v*S_v/S; % 3.193 (pg 160)
% Cnr = Cnr_wf + Cnr_v; % 3.192 (pg 160)

%% Trim Estimates & Assumptions
det = CLa*CMde-CMa*CLde;
a1 = ((CL1-CL0)*CMde+(CM0)*CLde)/(det)*180/pi;
de = (-CLa*(CM0)-CMa*(CL1-CL0))/det*180/pi;
CD1 = -0.0039*a1+0.0053; % nolinear fit from wind tunnel data
CDa = -0.0039*180/pi;
CDde = 0.00058*180/pi; % 1/rad,this is a linear fit
Ctxu = 0; % Assume thrust dose not change with speed
CMt1 = 0; % Assume thrust does not creat any pitching moment
CMta = 0;
CMtu = 0;

%% Thrust
T=1*CD1*qbar*S; % thrust at cruise condition is assumed to equal to drag
Ctx1=T/(qbar*S);

%% Dimesional Stability Derivatives
Xu = -qbar*S*(CDu+2*CD1)/(m*U1); %need to check
Xtu = qbar*S*(Ctxu+2*Ctx1)/(m*U1); % 
Xa = -qbar*S*(CDa-CL1)/m;
Xde = -qbar*S*CDde/m;% 
Zu = -qbar*S*(CLu+2*CL1)/(m*U1);% 
Za = -qbar*S*(CLa+CD1)/m; %
Zadot = -qbar*S*cbar*CLadot/(2*m*U1);% 
Zq = -qbar*S*cbar*CLq/(2*m*U1);% 
Zde = -qbar*S*CLde/m; %
Mu = qbar*S*cbar*(CMu+2*CM1)/(IYY*U1); % 
Mtu = qbar*S*cbar*(CMtu+CMt1)/(IYY*U1); % Cmtu,Cmt1
Ma = qbar*S*cbar*CMa/IYY; %
Mta = qbar*S*cbar*CMta/IYY; % Cmta
Madot = qbar*S*cbar^2*CMadot/(2*IYY*U1);% 
Mq = qbar*S*cbar^2*CMq/(2*IYY*U1);% 
Mde = qbar*S*cbar*CMde/IYY;% 

%% State Space Model
% E xdot = A x + B u
% states:
% -------
% alpha -> angle of attack
% q     -> pitch rate
% vt    -> true airspeed
% theta -> pitch angle


Abar = [   Za     U1+Zq   Zu-Xtu*sin(a1+at)   -g*sin(gammae)
       Ma+Mta    Mq         Mu+Mtu               0
        Xa       0     Xu+Xtu*cos(a1+at)   -g*cos(gammae)
        0        1            0                  0        ];

Bbar = [ Zde
      Mde
      Xde
       0  ];

Ebar = [U1-Zadot   0   0   0
     -Madot     1   0   0
       0        0   1   0
       0        0   0   1];



