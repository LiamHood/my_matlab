clc; clear; close all;

%% Group 1
G1_sun = xlsread('Group 1 Part 2.xlsx',1);

%sun trial 1
t_sun_1 = G1_sun(:,1); %time for trial 1 (s)
p_al_sun = G1_sun(:,2); %polished aluminum temp (C)
b_al_s_sun = G1_sun(:,3); %small black aluminum temp
b_b_sun = G1_sun(:,4); %black brass temp

%sun trial 2
t_sun_2 = G1_sun(:,7); %time for trial 2 (2)
b_al_l_sun = G1_sun(:,8); %large black aluminum temp
w_al_sun = G1_sun(:,9); %while aluminum temp

G1_lamp = xlsread('Group 1 Part 2.xlsx',2);

%lamp trial 1
t_lamp_1 = G1_lamp(:,1); %time for trial 1 (s)
p_al_lamp = G1_lamp(:,2); %polished aluminum temp (C)
b_al_s_lamp = G1_lamp(:,3); %small black aluminum temp
b_b_lamp = G1_lamp(:,4); %black brass temp

%lamp trial 2
t_lamp_2 = G1_lamp(:,7); %time for trial 2 (2)
b_al_l_lamp = G1_lamp(:,8); %large black aluminum temp
w_al_lamp = G1_lamp(:,9); %while aluminum temp

E = ones(length(t_sun_1),1);

%finding slopes (sun)
%p_al: linear plot from t = 480:720
pf_p_al_sun = polyfit(t_sun_1(8:12),p_al_sun(8:12),1);
dTdt_p_al_sun = pf_p_al_sun(1); %change in deg C per sec

%b_al_s: linear plot from t = 240:420
pf_b_al_s_sun = polyfit(t_sun_1(4:7),b_al_s_sun(4:7),1);
dTdt_b_al_s_sun = pf_b_al_s_sun(1); %change in deg C per sec

%b_b: linear plot from t = 360:480
pf_b_b_sun = polyfit(t_sun_1(6:8),b_b_sun(6:8),1);
dTdt_b_b_sun = pf_b_b_sun(1); %change in deg C per sec

%b_al_l: linear plot from t = 540:720
pf_b_al_l_sun = polyfit(t_sun_1(9:12),b_al_l_sun(9:12),1);
dTdt_b_al_l_sun = pf_b_al_l_sun(1); %change in deg C per sec

%w_al: linear plot from t = 900:1140
pf_w_al_sun = polyfit(t_sun_1(15:19),w_al_sun(15:19),1);
dTdt_w_al_sun = pf_w_al_sun(1); %change in deg C per sec

%finding slopes (lamp)
%p_al: linear plot from t = 600:720
pf_p_al_lamp = polyfit(t_sun_1(10:12),p_al_lamp(10:12),1);
dTdt_p_al_lamp = pf_p_al_lamp(1); %change in deg C per sec

%b_al_s: linear plot from t = 180:300
pf_b_al_s_lamp = polyfit(t_sun_1(3:5),b_al_s_lamp(3:5),1);
dTdt_b_al_s_lamp = pf_b_al_s_lamp(1); %change in deg C per sec

%b_b: linear plot from t = 300:360
pf_b_b_lamp = polyfit(t_sun_1(5:6),b_b_lamp(5:6),1);
dTdt_b_b_lamp = pf_b_b_lamp(1); %change in deg C per sec

%b_al_l: linear plot from t = 240:360
pf_b_al_l_lamp = polyfit(t_sun_1(9:12),b_al_l_lamp(9:12),1);
dTdt_b_al_l_lamp = pf_b_al_l_lamp(1); %change in deg C per sec

%w_al: linear plot from t = 360:480
pf_w_al_lamp = polyfit(t_sun_1(6:8),w_al_lamp(6:8),1);
dTdt_w_al_lamp = pf_w_al_lamp(1); %change in deg C per sec

m = [34.14 29.46 91.4 52.62 34.14]; %(g)order: p, b_al_s, b_b, b_al_l ,w
m = m/1000; %mass in kg
c_al = 921; %J/kgC
c_brass = 920; %J/kgC

I_p_al_sun = (m(1)*c_al*dTdt_p_al_sun)/(0.0381*0.01905);
I_b_al_s_sun = (m(2)*c_al*dTdt_b_al_s_sun)/(0.0381*0.01905);
I_b_b_sun = (m(3)*c_brass*dTdt_b_b_sun)/(0.0381*0.01905);
I_b_al_l_sun = (m(4)*c_al*dTdt_b_al_l_sun)/(0.0381*0.0254);
I_w_al_sun = (m(5)*c_al*dTdt_w_al_sun)/(0.0381*0.01905);

I_p_al_lamp = (m(1)*c_al*dTdt_p_al_lamp)/(0.0381*0.01905);
I_b_al_s_lamp = (m(2)*c_al*dTdt_b_al_s_lamp)/(0.0381*0.01905);
I_b_b_lamp = (m(3)*c_brass*dTdt_b_b_lamp)/(0.0381*0.01905);
I_b_al_l_lamp = (m(4)*c_al*dTdt_b_al_l_lamp)/(0.0381*0.0254);
I_w_al_lamp = (m(5)*c_al*dTdt_w_al_lamp)/(0.0381*0.01905);


figure
errorbar(t_sun_1,p_al_sun,E,'*')
hold on
errorbar(t_sun_1,b_al_s_sun,E,'*')
errorbar(t_sun_1,b_al_l_sun,E,'*')
errorbar(t_sun_1,w_al_sun,E,'*')
errorbar(t_sun_1,b_b_sun,E,'*')
grid on
xlabel('Time (s)')
ylabel('Temp (C)')
title('Cylinder Temp vs. Time in Sun (Group 1)')
legend('Polished Al','Black Al (.75 in)','Black Al (1 in)',...
    'White Al','Black Brass','location','SouthEast')


E = ones(length(t_lamp_1),1);

figure
errorbar(t_lamp_1,p_al_lamp,E,'*')
hold on
errorbar(t_lamp_1,b_al_s_lamp,E,'*')
errorbar(t_lamp_1,b_al_l_lamp,E,'*')
errorbar(t_lamp_1,w_al_lamp,E,'*')
errorbar(t_lamp_1,b_b_lamp,E,'*')
grid on
xlabel('Time (s)')
ylabel('Temp (C)')
title('Cylinder Temp vs. Time Under Lamp (Group 1)')
legend('Polished Al','Black Al (.75 in)','Black Al (1 in)',...
    'White Al','Black Brass','location','SouthEast')
%% Group 2
G2_1 = xlsread('group2_cylinder_1.xlsx');
G2_2 = xlsread('group2_cylinder_2.csv');

t_G2_1 = G2_1(:,1);
b_al_s_lamp = G2_1(:,2);
w_al_lamp = G2_1(:,3);
p_al_lamp = G2_1(:,4);

t_G2_2 = G2_2(:,1);
b_b_lamp = G2_2(:,2);
b_al_l_lamp = G2_2(:,3);

%already calculated by group 2
dTdt_b_al_s_lamp(2) = 0.044;
dTdt_w_al_lamp(2) = 0.036;
dTdt_p_al_lamp(2) = 0.014;
dTdt_b_b_lamp(2) = 0.0306;
dTdt_b_al_l_lamp(2) = 0.05;

I_p_al_lamp(2) = (m(1)*c_al*dTdt_p_al_lamp(2))/(0.0381*0.01905);
I_b_al_s_lamp(2) = (m(2)*c_al*dTdt_b_al_s_lamp(2))/(0.0381*0.01905);
I_b_b_lamp(2) = (m(3)*c_brass*dTdt_b_b_lamp(2))/(0.0381*0.01905);
I_b_al_l_lamp(2) = (m(4)*c_al*dTdt_b_al_l_lamp(2))/(0.0381*0.0254);
I_w_al_lamp(2) = (m(5)*c_al*dTdt_w_al_lamp(2))/(0.0381*0.01905);

E_1 = ones(length(t_G2_1),1)*.5;
E_2 = ones(length(t_G2_2),1)*.5;

figure
errorbar(t_G2_1,p_al_lamp,E_1)
hold on
errorbar(t_G2_1,b_al_s_lamp,E_1)
errorbar(t_G2_2,b_al_l_lamp,E_2)
errorbar(t_G2_1,w_al_lamp,E_1)
errorbar(t_G2_2,b_b_lamp,E_2)
grid on
xlabel('Time (s)')
ylabel('Temp (C)')
title('Cylinder Temp vs. Time Under Lamp (Group 2)')
legend('Polished Al','Black Al (.75 in)','Black Al (1 in)',...
    'White Al','Black Brass','location','SouthEast')

%% Group 3
G3_1 = xlsread('3_aluminium_all_finishes_heatup.csv');
G3_2 = xlsread('Brass_and_Large_Al_heatup_and_cooldown.csv');

t3_1 = G3_1(:,1);
b_al_s_lamp = G3_1(:,2);
b_al_l_lamp = G3_1(:,3);
w_al_lamp = G3_1(:,4);

t3_2 = G3_2(:,1);
p_al_lamp = G3_2(:,26);
b_b_lamp = G3_2(:,27);

E_1 = ones(length(t3_1),1)*.5;
E_2 = ones(length(t3_2),1)*.5;

figure
errorbar(t3_2,p_al_lamp,E_2)
hold on
errorbar(t3_1,b_al_s_lamp,E_1)
errorbar(t3_1,b_al_l_lamp,E_1)
errorbar(t3_1,w_al_lamp,E_1)
errorbar(t3_2,b_b_lamp,E_2)
grid on
xlabel('Time (s)')
ylabel('Temp (C)')
title('Cylinder Temp vs. Time Under Lamp (Group 3)')
legend('Polished Al','Black Al (.75 in)','Black Al (1 in)',...
    'White Al','Black Brass','location','NorthWest')

%finding slopes (lamp)
%p_al: linear plot from t = 600:720
pf_p_al_lamp = polyfit(t3_2(500:1000),p_al_lamp(500:1000),1);
dTdt_p_al_lamp(3) = pf_p_al_lamp(1); %change in deg C per sec

%b_al_s: linear plot from t = 180:300
pf_b_al_s_lamp = polyfit(t3_1(300:500),b_al_s_lamp(300:500),1);
dTdt_b_al_s_lamp(3) = pf_b_al_s_lamp(1); %change in deg C per sec

%b_b: linear plot from t = 300:360
pf_b_b_lamp = polyfit(t3_2(500:1000),b_b_lamp(500:1000),1);
dTdt_b_b_lamp(3) = pf_b_b_lamp(1); %change in deg C per sec

%b_al_l: linear plot from t = 240:360
pf_b_al_l_lamp = polyfit(t3_1(300:500),b_al_l_lamp(300:500),1);
dTdt_b_al_l_lamp(3) = pf_b_al_l_lamp(1); %change in deg C per sec

%w_al: linear plot from t = 360:480
pf_w_al_lamp = polyfit(t3_1(300:500),w_al_lamp(300:500),1);
dTdt_w_al_lamp(3) = pf_w_al_lamp(1); %change in deg C per sec

I_b_al_s_lamp(3) = (m(2)*c_al*dTdt_b_al_s_lamp(3))/(0.0381*0.01905);
I_b_b_lamp(3) = (m(3)*c_brass*dTdt_b_b_lamp(3))/(0.0381*0.01905);
I_b_al_l_lamp(3) = (m(4)*c_al*dTdt_b_al_l_lamp(3))/(0.0381*0.0254);
I_w_al_lamp(3) = (m(5)*c_al*dTdt_w_al_lamp(3))/(0.0381*0.01905);

%% Group 4
G4_1 = xlsread('Cylinders Group 4.xlsx',1);
G4_2 = xlsread('Cylinders Group 4.xlsx',2);

t4_1 = G4_1(:,1);
b_al_l_sun = G4_1(:,2);
w_al_sun = G4_1(:,3);
p_al_sun = G4_1(:,4);
b_b_sun = G4_1(:,6);
b_al_s_sun = G4_1(:,7);

t4_2 = G4_2(:,1);
b_al_l_lamp = G4_2(:,2);
w_al_lamp = G4_2(:,3);
p_al_lamp = G4_2(:,4);
b_b_lamp = G4_2(:,6);
b_al_s_lamp = G4_2(:,7);

E_1 = ones(length(t4_1),1)*.5;
E_2 = ones(length(t4_2),1)*.5;

figure
errorbar(t4_1,p_al_sun,E_1,'*')
hold on
errorbar(t4_1,b_al_s_sun,E_1,'*')
errorbar(t4_1,b_al_l_sun,E_1,'*')
errorbar(t4_1,w_al_sun,E_1,'*')
errorbar(t4_1,b_b_sun,E_1,'*')
grid on
xlabel('Time (s)')
ylabel('Temp (C)')
title('Cylinder Temp vs. Time in Sun (Group 4)')
legend('Polished Al','Black Al (.75 in)','Black Al (1 in)',...
    'White Al','Black Brass','location','NorthWest')

figure
errorbar(t4_2,p_al_lamp,E_2,'*')
hold on
errorbar(t4_2,b_al_s_lamp,E_2,'*')
errorbar(t4_2,b_al_l_lamp,E_2,'*')
errorbar(t4_2,w_al_lamp,E_2,'*')
errorbar(t4_2,b_b_lamp,E_2,'*')
grid on
xlabel('Time (s)')
ylabel('Temp (C)')
title('Cylinder Temp vs. Time Under Lamp (Group 4)')
legend('Polished Al','Black Al (.75 in)','Black Al (1 in)',...
    'White Al','Black Brass','location','NorthWest')


%finding slopes (sun)
%p_al: linear plot from t = 480:720
pf_p_al_sun = polyfit(t4_1(10:70),p_al_sun(10:70),1);
dTdt_p_al_sun(2) = pf_p_al_sun(1); %change in deg C per sec

%b_al_s: linear plot from t = 240:420
pf_b_al_s_sun = polyfit(t4_1(10:70),b_al_s_sun(10:70),1);
dTdt_b_al_s_sun(2) = pf_b_al_s_sun(1); %change in deg C per sec

%b_b: linear plot from t = 360:480
pf_b_b_sun = polyfit(t4_1(10:70),b_b_sun(10:70),1);
dTdt_b_b_sun(2) = pf_b_b_sun(1); %change in deg C per sec

%b_al_l: linear plot from t = 540:720
pf_b_al_l_sun = polyfit(t4_1(10:70),b_al_l_sun(10:70),1);
dTdt_b_al_l_sun(2) = pf_b_al_l_sun(1); %change in deg C per sec

%w_al: linear plot from t = 900:1140
pf_w_al_sun = polyfit(t4_1(10:70),w_al_sun(10:70),1);
dTdt_w_al_sun(2) = pf_w_al_sun(1); %change in deg C per sec

%finding slopes (lamp)
%p_al: linear plot from t = 600:720
pf_p_al_lamp = polyfit(t4_2(10:70),p_al_lamp(10:70),1);
dTdt_p_al_lamp(4) = pf_p_al_lamp(1); %change in deg C per sec

%b_al_s: linear plot from t = 180:300
pf_b_al_s_lamp = polyfit(t4_2(10:70),b_al_s_lamp(10:70),1);
dTdt_b_al_s_lamp(4) = pf_b_al_s_lamp(1); %change in deg C per sec

%b_b: linear plot from t = 300:360
pf_b_b_lamp = polyfit(t4_2(10:70),b_b_lamp(10:70),1);
dTdt_b_b_lamp(4) = pf_b_b_lamp(1); %change in deg C per sec

%b_al_l: linear plot from t = 240:360
pf_b_al_l_lamp = polyfit(t4_2(10:70),b_al_l_lamp(10:70),1);
dTdt_b_al_l_lamp(4) = pf_b_al_l_lamp(1); %change in deg C per sec

%w_al: linear plot from t = 360:480
pf_w_al_lamp = polyfit(t4_2(10:70),w_al_lamp(10:70),1);
dTdt_w_al_lamp(4) = pf_w_al_lamp(1); %change in deg C per sec

I_p_al_sun(2) = (m(1)*c_al*dTdt_p_al_sun(2))/(0.0381*0.01905);
I_b_al_s_sun(2) = (m(2)*c_al*dTdt_b_al_s_sun(2))/(0.0381*0.01905);
I_b_b_sun(2) = (m(3)*c_brass*dTdt_b_b_sun(2))/(0.0381*0.01905);
I_b_al_l_sun(2) = (m(4)*c_al*dTdt_b_al_l_sun(2))/(0.0381*0.0254);
I_w_al_sun(2) = (m(5)*c_al*dTdt_w_al_sun(2))/(0.0381*0.01905);

I_p_al_lamp(4) = (m(1)*c_al*dTdt_p_al_lamp(4))/(0.0381*0.01905);
I_b_al_s_lamp(4) = (m(2)*c_al*dTdt_b_al_s_lamp(4))/(0.0381*0.01905);
I_b_b_lamp(4) = (m(3)*c_brass*dTdt_b_b_lamp(4))/(0.0381*0.01905);
I_b_al_l_lamp(4) = (m(4)*c_al*dTdt_b_al_l_lamp(4))/(0.0381*0.0254);
I_w_al_lamp(4) = (m(5)*c_al*dTdt_w_al_lamp(4))/(0.0381*0.01905);

mean_I_p_al_lamp = mean(I_p_al_lamp);
mean_I_b_al_l_lamp = mean(I_b_al_l_lamp);
mean_I_b_al_s_lamp = mean(I_b_al_s_lamp);
mean_I_b_b_lamp = mean(I_b_b_lamp);
mean_I_w_al_lamp = mean(I_w_al_lamp);

mean_I_p_al_sun = mean(I_p_al_sun);
mean_I_b_al_l_sun = mean(I_b_al_l_sun);
mean_I_b_al_s_sun = mean(I_b_al_s_sun);
mean_I_b_b_sun = mean(I_b_b_sun);
mean_I_w_al_sun = mean(I_w_al_sun);

std_I_p_al_lamp = std(I_p_al_lamp);
std_I_b_al_l_lamp = std(I_b_al_l_lamp);
std_I_b_al_s_lamp = std(I_b_al_s_lamp);
std_I_b_b_lamp = std(I_b_b_lamp);
std_I_w_al_lamp = std(I_w_al_lamp);

std_I_p_al_sun = std(I_p_al_sun);
std_I_b_al_l_sun = std(I_b_al_l_sun);
std_I_b_al_s_sun = std(I_b_al_s_sun);
std_I_b_b_sun = std(I_b_b_sun);
std_I_w_al_sun = std(I_w_al_sun);

