%% Orbits Project 2 - Jason Dove - Liam Hood
% Initial Vectorsc
clear all; close all; clc ;
TimeStart = [2019,11,21,12,0,0] ;
JDStart = juliandate(TimeStart) ;
global muearth REarth mu
muearth = 398600 ; % km^3/s^2
mu = 398600 ; % km^3/s^2
REarth = 6378 ; % km - Radius of Earth

%% Satellite: FLYING LAPTOP           
% Catalog Number: 42831
Epochtime1 = 19308.91755737 ;
Inc1 =  97.545900 ; %deg
RAAN1 =  195.955300 ; %deg
Ecc1 = 0.001324 ; 
ARGOPER1 = 322.069700 ; %deg
MeanAnom1 = 37.958900 ; %deg
MeanMot1 = 14.910514 ; %rev/day
% Period of rev: 5795 s/rev
% Semi-major axis: 6972790 meters
% Semi-minor axis: 6972784 meters

State1 = TLE2RV(TimeStart,Epochtime1,muearth,MeanMot1,MeanAnom1,Ecc1,Inc1,RAAN1,ARGOPER1) ;

%% Satellite: VIASAT-1                
% Catalog Number: 37843
Epochtime2 = 19309.54287057 ;
Inc2 = 0.018600 ; %deg
RAAN2 =  265.075300 ; %deg
Ecc2 = 0.000226 ;
ARGOPER2 =  311.649000 ; % deg
MeanAnom2 = 268.084300 ; %deg
MeanMot2 = 1.002724 ; % rev/day
% Period of rev: 86165 s/rev
% Semi-major axis: 42164539 meters
% Semi-minor axis: 42164538 meters

State2 = TLE2RV(TimeStart,Epochtime2,muearth,MeanMot2,MeanAnom2,Ecc2,Inc2,RAAN2,ARGOPER2) ;

%% STARLINK-31             
% 1 44235U 19029A   19339.50001157  .00034642  00000-0  19953-2 0  9991
% 2 44235  53.0039  14.8578 0001329 117.2061 179.6033 15.12119937 30726
% Satellite: STARLINK-31             
% Catalog Number: 44235
% Epoch time: 19339.50001157
% Inclination: 53.003900 deg
% RA of ascending node: 14.857800 deg
% Eccentricity: 0.000133
% Arg of perigee: 117.206100 deg
% Mean anomaly: 179.603300 deg
% Mean motion: 15.121199 rev/day
% Period of rev: 5714 s/rev
% Semi-major axis: 6907871 meters
% Semi-minor axis: 6907870 meters
Epochtime3 = 19339.50001157 ;
Inc3 = 53.003900 ; %deg
RAAN3 = 14.857800; %.784600 deg
Ecc3 = 0.000133 ;
ARGOPER3 = 117.206100 ; %deg
MeanAnom3 = 179.603300 ; %deg
MeanMot3 =  15.121199 ; %rev/day
% Period of rev: 5714 s/rev
% Semi-major axis: 7101002 meters
% Semi-minor axis: 7101002 meters

State3 = TLE2RV(TimeStart,Epochtime3,muearth,MeanMot3,MeanAnom3,Ecc3,Inc3,RAAN3,ARGOPER3) ;
% Radius of Apogee and Perigee Functions
RP = @(h,ecc) h.^2./(muearth*(1+ecc)) ;
RA = @(h,ecc) h.^2./(muearth*(1-ecc)) ;
%% FL Run Over Certain Period No Correction
options = odeset('Events',@MyEvent,'RelTol' , 1e-8 , 'AbsTol' , 1e-8) ;
% Flying Laptop
TimePassed = 60*60*24*20 ; %2 Days in Seconds
NBody1 = [1;0] ; Drag1 = 1 ; JVec1 = ones(1,6) ; SRP1 = 1 ; Area1 = 1 ; Mass1 = 120 ;
[tnew1 , StateNew1,~] = ode45(@GravODE,[0,TimePassed],[State1;State1],options,muearth,NBody1,Drag1,JVec1,SRP1,JDStart,Area1,Mass1) ;
% [tnew1,StateNew1,Burns1,BurnTime1] = CorrectBurn(State1,muearth,TimePassed,JDStart,Area1,Mass1,Inf) ;
[h1,inc1,ecc1,RAAN1,ARGOPER1,~,~] = orbitvars(StateNew1(:,1:3)',StateNew1(:,4:6)',muearth) ;
RP1 = RP(h1,ecc1) ;
RA1 = RA(h1,ecc1) ;


%% VIASAT-1 Run Over Certain Period Before Correction
options = odeset('Events',@MyEvent,'RelTol' , 1e-8 , 'AbsTol' , 1e-8) ;
% VIASAT-1
TimePassed = 60*60*24*365*5; %2 Days in Seconds
HydArea = @(X,Y,Z) pi*(3/4/pi*(X*Y*Z))^(2/3) ; % Hydraulic Area of Object
Area2 = HydArea(3.409,8.733,3.061); Mass2 = 6740 ;
NBody2 = [1;0] ; Drag2 = 1 ; JVec2 = ones(1,6) ; SRP2 = 1 ;  %Need Changed Area
% [tnew2 , StateNew2,~] = ode45(@GravODE,[0,TimePassed],[State2;State2],options,muearth,NBody2,Drag2,JVec2,SRP2,JDStart,Area2,Mass2) ;
[tnew2,StateNew2,Burns2,BurnTime2] = CorrectBurn(State2,muearth,TimePassed,JDStart,Area2,Mass2,Inf) ;
[h2,inc2,ecc2,RAAN2,ARGOPER2,~,~] = orbitvars(StateNew2(:,1:3)',StateNew2(:,4:6)',muearth) ;
RP2 = RP(h2,ecc2) ;
RA2  = RA(h2,ecc2) ;

%% Starlink Sat-One Correction
TimePassed = 60*60*24*10; %10 Days in Seconds
HydArea = @(X,Y,Z) pi*(3/4/pi*(X*Y*Z))^(2/3) ; % Hydraulic Area of Object
Area4 = .07^2+8*2 ; Mass4 = 227 ;
NBody4 = [1;0] ; Drag4 = 1 ; JVec4 = ones(1,6) ; SRP4 = 1 ;  %Need Changed Area
% [tnew4 , StateNew4,~] = ode45(@GravODE,[0,TimePassed],[State3;State3],options,muearth,NBody4,Drag4,JVec4,SRP4,JDStart,Area4,Mass4) ;
[tnew4,StateNew4,Burns4,BurnTime4] = CorrectBurn(State3,muearth,TimePassed,JDStart,Area4,Mass4,Inf) ;
[h4,inc4,ecc4,RAAN4,ARGOPER4,~,~] = orbitvars(StateNew4(:,1:3)',StateNew4(:,4:6)',muearth) ;
RP4 = RP(h4,ecc4) ;
RA4  = RA(h4,ecc4) ;
%% Continuous Correction of Viasat
% Viasat
TimePassed = 60*60*24*365*5; %2 Days in Seconds
Area3 = HydArea(3.409,8.733,3.061) ; Mass3 = 6740 ;
NBody3 = [0;0] ; Drag3 = 0 ; JVec3 = ones(1,6)*0 ; SRP3 = 0 ;  %Need Changed Area
[tnew3NoP , StateNew3NoP,~] = ode45(@GravODE,[0,TimePassed],[State2;State2],options,muearth,NBody3,Drag3,JVec3,SRP3,JDStart,Area3,Mass3) ;
[tnew3,StateNew3,Burns3,BurnTime3] = CorrectBurn(State2,muearth,TimePassed,JDStart,Area3,Mass3,1000) ;
[h3,inc3,ecc3,RAAN3,ARGOPER3,~,~] = orbitvars(StateNew3(:,1:3)',StateNew3(:,4:6)',muearth) ;
RP3 = RP(h3,ecc3) ;
RA3  = RA(h3,ecc3) ;
%% Extra
% Viasat
TimePassed = 60*60*24*365*5; %2 Days in Seconds
Area5 = HydArea(3.409,8.733,3.061) ; Mass5 = 6740 ;
NBody5 = [1;0] ; Drag5 = [1,1] ; JVec5 = ones(1,6) ; SRP5 = [1,1] ;  %Need Changed Area
[tnew5 , StateNew5,~] = ode45(@GravODE,[0,TimePassed],[State2;State2],options,muearth,NBody5,Drag5,JVec5,SRP5,JDStart,Area5,Mass5) ;
% NBody5 = [0;0] ; Drag5 = 0 ; JVec5 = ones(1,6)*0 ; SRP5 = 0 ;  %Need Changed Area
% [tnew5NoP , StateNew5NoP,~] = ode45(@GravODE,[0,TimePassed],[State2;State2],options,muearth,NBody5,Drag5,JVec5,SRP5,JDStart,Area5,Mass5) ;
% [tnew5,StateNew5,Burns5,BurnTime5] = CorrectBurn(State2,muearth,TimePassed,JDStart,Area5,Mass5,Inf) ;
[h5,inc5,ecc5,RAAN5,ARGOPER5,~,~] = orbitvars(StateNew5(:,1:3)',StateNew5(:,4:6)',muearth) ;
RP5 = RP(h5,ecc5) ;
RA5  = RA(h5,ecc5) ;
%% Plots For FL
figure 
plot3(StateNew1(:,1),StateNew1(:,2),StateNew1(:,3))
hold on
plot3(StateNew1(:,7),StateNew1(:,8),StateNew1(:,9))
axis equal
GlobeTrotter
grid on
legend({'Perturbed','Unperturbed'})
title('Flying Laptop: Orbit')

figure
plot(tnew1/60/60/24,RA1-REarth) 
hold on 
plot(tnew1/60/60/24,RP1-REarth) 
title('Flying Laptop: Perigee and Apogee Altitude vs. Time')
xlabel('Time [Days]')
ylabel('Altitude [km]')
legend({'Apogee','Perigee'})
Ylimits = get(gca,'YLim') ;
legend({'Perturbed Orbit','Unperturbed Orbit'})
% ylim([0,Ylimits(2)])
grid on

figure
subplot(4,1,1)
title('RAAN vs. Time')
plot(tnew1/60/60/24,RAAN1-RAAN1(1)) 
xlabel('Time [Days]')
ylabel('RAAN [degs]')
% ylim([0 360])
grid on
grid minor
title('Flying Laptop Element Changes')

subplot(4,1,2)
title('Inclination vs. Time')
plot(tnew1/60/60/24,inc1-inc1(1)) 
xlabel('Time [Days]')
ylabel('Inc [degs]')
% ylim([-180 180])
grid on
grid minor

subplot(4,1,3)
title('Argument of Perigee vs. Time')
plot(tnew1/60/60/24,ARGOPER1-ARGOPER1(1)) 
xlabel('Time [Days]')
ylabel('Omega [degs]')
% ylim([0 360])
grid on
grid minor

subplot(4,1,4)
title('Eccentricity vs. Time')
plot(tnew1/60/60/24,ecc1) 
xlabel('Time [Days]')
ylabel('Ecc')
% ylim([0 360])
grid on
grid minor
% suptitle('Flying Laptop Element Changes')
%% Dv
fprintf('Flying Laptop')
%% Plots for VIASAT-1
figure 
plot3(StateNew2(:,1),StateNew2(:,2),StateNew2(:,3))
axis equal
hold on
Pos2 = find(tnew2<BurnTime2(1),1,'last') ;
plot3(StateNew2(1:Pos2,7),StateNew2(1:Pos2,8),StateNew2(1:Pos2,9),'LineWidth',2)
GlobeTrotter
grid on

legend({'Perturbed','Unperturbed'})
title('VIASAT-1: Orbit w/ Corrective Burn')

figure
plot(tnew2/60/60/24,RA2-REarth) 
hold on 
plot(tnew2/60/60/24,RP2-REarth) 
title('VIASAT-1: Perigee and Apogee Altitude vs. Time w/ Corrective Burn')
xlabel('Time [Days]')
ylabel('Altitude [km]')
legend({'Apogee','Perigee'})
Ylimits = get(gca,'YLim') ;
% ylim([0,Ylimits(2)])
grid on

figure
subplot(4,1,1)
title('RAAN vs. Time')
plot(tnew2/60/60/24,RAAN2-RAAN2(1)) 
xlabel('Time [Days]')
ylabel('RAAN [degs]')
% ylim([0 360])
grid on
grid minor
title('VIASAT-1: Element Changes w/ Corrective Burn')

subplot(4,1,2)
title('Inclination vs. Time')
plot(tnew2/60/60/24,inc2-inc2(1)) 
xlabel('Time [Days]')
ylabel('Inc [degs]')
% ylim([-180 180])
grid on
grid minor

subplot(4,1,3)
title('Argument of Perigee vs. Time')
plot(tnew2/60/60/24,ARGOPER2-ARGOPER2(1)) 
xlabel('Time [Days]')
ylabel('Omega [degs]')
% ylim([0 360])
grid on
grid minor

subplot(4,1,4)
title('Eccentricity vs. Time')
plot(tnew2/60/60/24,ecc2) 
xlabel('Time [Days]')
ylabel('Ecc')
% ylim([0 360])
grid on
grid minor
% suptitle('VIASAT-1: Element Changes w/ Corrective Burn')

%% Dv
fprintf('VIASAT Single Correction')
fprintf('Correction dV is from a total of %s burns of with a summation of %fkm/s in dV.\n',num2str(length(Burns2)),sum(Burns2)) ;
%% Plots for Starlink-31
figure 
plot3(StateNew4(:,1),StateNew4(:,2),StateNew4(:,3))
axis equal
hold on
Pos4 = find(tnew4<BurnTime4(1),1,'last');
plot3(StateNew4(1:Pos4,7),StateNew4(1:Pos4,8),StateNew4(1:Pos4,9),'LineWidth',2)
GlobeTrotter
legend({'Perturbed','Unperturbed'})
grid on
title('STARLINK-31: Orbit w/ Corrective Burn')

figure
plot(tnew4/60/60/24,RA4-REarth) 
hold on 
plot(tnew4/60/60/24,RP4-REarth) 
title('STARLINK-31: Perigee and Apogee Altitude vs. Time w/ Corrective Burn')
xlabel('Time [Days]')
ylabel('Altitude [km]')
legend({'Apogee','Perigee'})
Ylimits = get(gca,'YLim') ;
% ylim([0,Ylimits(2)])
grid on

figure
subplot(4,1,1)
title('RAAN vs. Time')
plot(tnew4/60/60/24,RAAN4-RAAN4(1)) 
xlabel('Time [Days]')
ylabel('RAAN [degs]')
% ylim([0 360])
grid on
grid minor
title('STARLINK-31: Element Changes w/ Corrective Burn')

subplot(4,1,2)
title('Inclination vs. Time')
plot(tnew4/60/60/24,inc4-inc4(1)) 
xlabel('Time [Days]')
ylabel('Inc [degs]')
% ylim([-180 180])
grid on
grid minor

subplot(4,1,3)
title('Argument of Perigee vs. Time')
plot(tnew4/60/60/24,ARGOPER4-ARGOPER4(1)) 
xlabel('Time [Days]')
ylabel('Omega [degs]')
% ylim([0 360])
grid on
grid minor

subplot(4,1,4)
title('Eccentricity vs. Time')
plot(tnew4/60/60/24,ecc4) 
xlabel('Time [Days]')
ylabel('Ecc')
% ylim([0 360])
grid on
grid minor
% suptitle('STARLINK-31: Element Changes w/ Corrective Burn')
%% DV
fprintf('Starlink Single Correction')
fprintf('Correction dV is from a total of %s burns of with a summation of %fkm/s in dV.\n',num2str(length(Burns4)),sum(Burns4)) ;
%% Plots of Corrective VIASAT-1
figure 
plot3(StateNew3(:,1),StateNew3(:,2),StateNew3(:,3))
hold on
plot3(StateNew3NoP(:,7),StateNew3NoP(:,8),StateNew3NoP(:,9),'LineWidth',2)
axis equal
GlobeTrotter
grid on
legend({'Perturbed','Unperturbed'})
title({'Tolerance Bound Corrected Burns', 'ECI for VIASAT-1'})

figure
plot(tnew3/60/60/24,RA3-REarth) 
hold on 
plot(tnew3/60/60/24,RP3-REarth) 
title({'VIASAT-1: Perigee and Apogee Altitude vs. Time','Corrected Elements, Tolerance Bound'})
xlabel('Time [Days]')
ylabel('Altitude [km]')
legend({'Apogee','Perigee'})
Ylimits = get(gca,'YLim') ;
legend({'Apogee','Perigee'})
% ylim([0,Ylimits(2)])
grid on

figure
subplot(4,1,1)
title('RAAN vs. Time')
plot(tnew3/60/60/24,RAAN3-RAAN3(1)) 
xlabel('Time [Days]')
ylabel('RAAN [degs]')
% ylim([0 360])
grid on
grid minor
title({'VIASAT-1 Element Changes: Corrections' , 'Tolerance Bound'})

subplot(4,1,2)
title('Inclination vs. Time')
plot(tnew3/60/60/24,inc3-inc3(1)) 
xlabel('Time [Days]')
ylabel('Inc [degs]')
% ylim([-180 180])
grid on
grid minor

subplot(4,1,3)
title('Argument of Perigee vs. Time')
plot(tnew3/60/60/24,ARGOPER3-ARGOPER3(1)) 
xlabel('Time [Days]')
ylabel('Omega [degs]')
% ylim([0 360])
grid on
grid minor

subplot(4,1,4)
title('Eccentricity vs. Time')
plot(tnew3/60/60/24,ecc3) 
xlabel('Time [Days]')
ylabel('Ecc')
% ylim([0 360])
grid on
grid minor
% suptitle({'VIASAT-1 Element Changes: Corrections' , 'Tolerance Bound'})

%% DV
fprintf('VIASAT-1 Tolerance Bound Corrections')
fprintf('For a 1000 km Deviation Control: Correction dV is from a total of %s burns of with a summation of %fkm/s in dV.\n',num2str(length(Burns3)),sum(Burns3)) ;

%% Extra Sail  VIASAT-1 Plots
figure 
plot3(StateNew5(:,1),StateNew5(:,2),StateNew5(:,3))
hold on
plot3(StateNew5(:,7),StateNew5(:,8),StateNew5(:,9),'LineWidth',2)
axis equal
GlobeTrotter
grid on
legend({'Perturbed','Unperturbed'})
title({'Solar Sail Use', 'ECI for VIASAT-1'})

figure
plot(tnew5/60/60/24,RA5-REarth) 
hold on 
plot(tnew5/60/60/24,RP5-REarth) 
title({'VIASAT-1: Perigee and Apogee Altitude vs. Time','Solar Sail Use'})
xlabel('Time [Days]')
ylabel('Altitude [km]')
legend({'Apogee','Perigee'})
Ylimits = get(gca,'YLim') ;
legend({'Apogee','Perigee'})
% ylim([0,Ylimits(2)])
grid on

figure
subplot(4,1,1)
title('RAAN vs. Time')
plot(tnew5/60/60/24,RAAN5-RAAN5(1)) 
xlabel('Time [Days]')
ylabel('RAAN [degs]')
% ylim([0 560])
grid on
grid minor
title({'VIASAT-1 Element Changes: Solar Sail Use'})

subplot(4,1,2)
title('Inclination vs. Time')
plot(tnew5/60/60/24,inc5-inc5(1)) 
xlabel('Time [Days]')
ylabel('Inc [degs]')
% ylim([-180 180])
grid on
grid minor

subplot(4,1,3)
title('Argument of Perigee vs. Time')
plot(tnew5/60/60/24,ARGOPER5-ARGOPER5(1)) 
xlabel('Time [Days]')
ylabel('Omega [degs]')
% ylim([0 360])
grid on
grid minor

subplot(4,1,4)
title('Eccentricity vs. Time')
plot(tnew5/60/60/24,ecc5) 
xlabel('Time [Days]')
ylabel('Ecc')
% ylim([0 360])
grid on
grid minor
% suptitle({'VIASAT-1 Element Changes: Solar Sail Use'})

%% Corrective Burning
function [times,StateOut,Burns,BurnTime] = CorrectBurn(StateIn,muearth,Time,JDStart,Area2,Mass2,CorrectingDistance)
NBody2 = [1;0] ; Drag2 = 1 ; JVec2 = ones(1,6) ; SRP2 = 1 ; %Need Changed Area
hh = 0;
StateOut = [] ;
tOut = 0 ;
Burns = [] ;
BurnTime = [] ;
oo = 0;
while tOut(end) < Time
    options = odeset('Events',@MyDistEvent,'RelTol' , 1e-8 , 'AbsTol' , 1e-8) ;
    [tOutFirst , StateOutPart , ~ , ~ , Pos] = ode45(@GravODE,[tOut(end),Time],[StateIn;StateIn],options,muearth,NBody2,Drag2,JVec2,SRP2,JDStart,Area2,Mass2,CorrectingDistance) ;
    tOut = [tOut;tOutFirst] ;
    hh = hh + length(tOut) ;
    if ~isempty(Pos) == 1 || CorrectingDistance == Inf
        hek = 'ahhh' ;
        [ha1,inca1,ecca1,RAANa1,ARGOPER1,TA1,SMA] = orbitvars(StateOutPart(end,7:9)',StateOutPart(end,10:12)',muearth) ;
%         for ll = (1:359)
%             StateOutie(1:6,ll) = Orb2ECI(ha1,muearth,ecca1,TA1+ll,inca1,RAANa1,ARGOPER1) ;
%             Vec(ll) = dot(StateOutie(1:3,ll),StateOutPart(end,1:3)')/(norm(StateOutie(1:3,ll))*norm(StateOutPart(end,1:3)')) ;
%         end
%         [~ , Positioner] = min(Vec) ;
        HalfPeriod = pi*sqrt(SMA^3/muearth) ;
        for pp = (1:3590)
            StateOutie(1:6,pp) = Orb2ECI(ha1,muearth,ecca1,TA1+pp/10,inca1,RAANa1,ARGOPER1) ;
%             Vec(ll) = dot(StateOutie(1:3,ll),StateOutPart(end,1:3)')/(norm(StateOutie(1:3,ll))*norm(StateOutPart(end,1:3)')) ;
%             Timeso = linspace(HalfPeriod/100,HalfPeriod + HalfPeriod*3/4,1000) ;
            [V1a , V2a] = LambertoCalc(StateOutPart(end,1:3)',StateOutie(1:3,pp),HalfPeriod,muearth,1);
            if imag(V1a) == 0
                V1j(:,pp) = real(V1a) ;
                V2j(:,pp) = real(V2a) ;
            else
                V1j(:,pp) = [Inf;Inf;Inf] ;
                V2j(:,pp) = [Inf;Inf;Inf] ;
            end
        end
        for pp = 1:pp
            Burnse(pp) = norm(V1j(:,pp) - StateOutPart(end,4:6)') ;
        end
        [Burnses,posit] = min(abs(Burnse)) ;
    Burns(oo+1) = Burnses ;
    V1 = V1j(:,posit) ;
    V2 = V2j(:,posit) ;
    BurnTime(oo+1) = tOut(end) ;
    options = odeset('RelTol' , 1e-8 , 'AbsTol' , 1e-8) ;
    [tOutBurner , StateOutPartBurner] = ode45(@GravODE,[0,HalfPeriod],[StateOutPart(end,1:3)';V1;StateOutPart(end,7:12)'],options,muearth,[0;0],0,zeros(6,1),0,JDStart) ;
    tOut = [tOut;tOutBurner+tOut(end)] ;
    BurnTime(oo+2) = tOut(end) ;
    EnderState = StateOutie(1:6,posit)' ;
    StateOutPart = [StateOutPart;StateOutPartBurner(1:end-1,:);[EnderState,EnderState]] ;
    Burns(oo+2) = norm(StateOutPartBurner(end,4:6) - EnderState(4:6)) ;
    oo = oo + 2 ;
    end
    StateIn = StateOutPart(end,1:6)' ;
    StateOut = [StateOut;StateOutPart] ;
end
if isempty(Burns) == 1
    Burns = 0 ;
    BurnTime = 0 ;
end
times = tOut(2:end) ;
end
function  [value, isterminal, direction] = MyDistEvent(~, y,~,~,~,~,~,~,~,~,CorrectingDist)
value      =    norm(y(1:3) - y(7:9)) - CorrectingDist ;
isterminal =    1;   % Stop the integration
direction  =    0;
end
function [State] = Orb2ECI(h,mu,ecc,TA,inc,RAAN,ArgOfPerigee)
% Combines elements to perifocal and perifocal to ECI
PeriState = PERISTATE(h,mu,ecc,TA) ;
[State] = PeriToECI(PeriState,inc,RAAN,ArgOfPerigee) ;
function State = PERISTATE(h,mu,ecc,TA)
%input 4 values to get the state vector in a perifocal frame
State(1:3,1) = h^2/(mu*(1+ecc*cosd(TA)))*[cosd(TA);sind(TA);0] ; %position in perifocal

State(4:6,1) = mu/h*[-sind(TA);(ecc+cosd(TA));0] ; %Velocity in Perifocal
end
end
%% Perifocal to ECI
function [State] = PeriToECI(State,inc,RAAN,ArgOfPerigee)
%state must be input as a vector
    function [Mat] = R3(Angle) % generates a 3rd axis rotation matrix
        Mat = [cosd(Angle),-sind(Angle),0;sind(Angle),cosd(Angle),0;0,0,1] ;
    end
    function [Mat] = R1(Angle) % generates a 1st axis rotation matrix
        Mat = [1,0,0;0,cosd(Angle),-sind(Angle);0,sind(Angle),cosd(Angle)];
    end
State(1:3) = R3(RAAN)*R1(inc)*R3(ArgOfPerigee)*State(1:3) ; %transforms the position into perifocal
State(4:6) = R3(RAAN)*R1(inc)*R3(ArgOfPerigee)*State(4:6) ; %transforms the velocity into perifocal
end

%% ODE Propogation
function [ dalldt ] = GravODE(t,state,mu,NBody,Drag,JVec,SRPA,jd,Area,Mass,CorrectingDist)
% Find dAcceleration 
if length(SRPA)== 2
    SRPONOFF = 1 ;
else
    SRPONOFF = 0 ;
end
SRP = SRPA(1) ;
if nargin == 8 
    Area = Inf ;
    Mass = Inf ;
    CorrectingDist = Inf ;
end
if nargin < 11
    CorrectingDist = Inf;
end
r = state(1:3) ;
v = state(4:6) ;
rNoP = state(7:9) ;
vNoP = state(10:12) ;
AccNoP = -mu*rNoP/(norm(rNoP))^3 ;
if length(NBody) >= 1
    p = NBod(r,NBody,t,jd) ;
else
    p = 0 ;
end
if Drag == 1
    p = p +  Dragger([r;v],2.2,Area,Mass) ;
end
if SRP == 1 
    p = p + SRPacceleration( r , v , Area , Mass , jd , t , SRPONOFF ) ;
end
if length(JVec) ~= 1
    p = p + J26(state(1:6),JVec,mu) ;
end
Acc = -mu*r/(norm(r))^3 + p ; % Acc is in km 
% Combine ( first Part is the velocity change is acc )
dalldt = [v;Acc;vNoP;AccNoP] ;
end

%% Code For J2-J6 Perturbations
function p = J26(State,Js,mu)
r = State(1:3) ;
v = State(4:6) ;
rnorm = norm(r) ;
REarth = 6378 ; % km - Radius of Earth
%Perturbation
p = 0 ;
if Js(1) == 1 % J2 Addition
    J2 =  1.08262668355E-3 ;
    apil = -3*J2*mu*REarth^2*r(1)/(2*rnorm^5)*(1-5*r(3)^2/rnorm^2);
    apjl = -3*J2*mu*REarth^2*r(2)/(2*rnorm^5)*(1-5*r(3)^2/rnorm^2);
    apkl = -3*J2*mu*REarth^2*r(3)/(2*rnorm^5)*(3-5*r(3)^2/rnorm^2);
    p = p + [apil;apjl;apkl] ;
end
if Js(2) == 1 % J3 Addition
    J3 = -2.53265648533E-6 ;
    api2 = -5*J3*mu*REarth^3*r(1)/(2*rnorm^7)*(3*r(3)-7*r(3)^3/rnorm^2);
    apj2 = -5*J3*mu*REarth^3*r(2)/(2*rnorm^7)*(3*r(3)-7*r(3)^3/rnorm^2);
    apk2 = -5*J3*mu*REarth^3/(2*rnorm^7)*(6*r(3)^2-7*r(3)^4/rnorm^2-3/5*rnorm^2);
    p = p + [api2;apj2;apk2] ;
end
if Js(3) == 1 % J4 Addition
    J4 = -1.61962159137E-6 ;
    api3 = 15*J4*mu*REarth^4*r(1)/(8*rnorm^7)*(1-14*r(3)^2/rnorm^2+21*r(3)^4/rnorm^4);
    apj3 = 15*J4*mu*REarth^4*r(2)/(8*rnorm^7)*(1-14*r(3)^2/rnorm^2+21*r(3)^4/rnorm^4);
    apk3 = 15*J4*mu*REarth^4*r(3)/(8*rnorm^7)*(5-70*r(3)^2/(3*rnorm^2)+21*r(3)^4/rnorm^4);
    p = p + [api3;apj3;apk3] ;
end
if Js(4) == 1 % J5 Addition
    J5 = -2.27296082869E-7 ;
    api4 = 3*J5*mu*REarth^5*r(1)*r(3)/(8*rnorm^9)*(35-210*r(3)^2/rnorm^2+231*r(3)^4/rnorm^4);
    apj4 = 3*J5*mu*REarth^5*r(2)*r(3)/(8*rnorm^9)*(35-210*r(3)^2/rnorm^2+231*r(3)^4/rnorm^4);
    apk4 = 3*J5*mu*REarth^5*r(3)^2/(8*rnorm^9)*(105-315*r(3)^2/(rnorm^2)+231*r(3)^4/rnorm^4)-15*J5*mu*REarth^5/(8*rnorm^7);
    p = p + [api4;apj4;apk4] ;
end
if Js(5) == 1 % J6 Addition
    J6 =  5.40681239107E-7 ;
    api5 = -J6*mu*REarth^6*r(1)/(16*rnorm^9)*(35-945*r(3)^2/rnorm^2+3465*r(3)^4/rnorm^4-3003*r(3)^6/rnorm^6);
    apj5 = -J6*mu*REarth^6*r(2)/(16*rnorm^9)*(35-945*r(3)^2/rnorm^2+3465*r(3)^4/rnorm^4-3003*r(3)^6/rnorm^6);
    apk5 = -J6*mu*REarth^6*r(3)/(16*rnorm^9)*(245-2205*r(3)^2/rnorm^2+4851*r(3)^4/rnorm^4-3003*r(3)^6/rnorm^6);
    p = p + [api5;apj5;apk5] ;
end
end
%% Liam SRP Codes
function [ light , us ] = ValladoShadow( rsc , JDo , t )
% [ light , rsun ] = ValladoShadow( rsc , JDo , t )
   radEarth = 6378 ;
   alphau = .264121687*(pi/180) ;
   alphap = .269007205*(pi/180) ;   
   
   JD = JDo + t/( 24*60*60 ) ;
   n = JD - 2451545.0 ;
   rbavg = 149597870.691 ;
   M = 357.529 + 0.98560023*n ;
   L = 280.459 + 0.98564736*n ;
   lam = L + 1.915*sind( M ) + 0.0200*sind( 2*M ) ;
   eps = 23.439 - ( 3.56e-7 ) * n ;
   us = [ cosd( lam ) ; sind( lam )*cosd( eps ) ; sind( lam )*sin( eps ) ] ;
   rbmag = ( 1.00014 - 0.01671*cosd( M ) - 0.000140*cos( 2*M ) )*rbavg ;
   rsun = rbmag*us ;
   rsunmag = rbmag ;
   rscmag = norm( rsc ) ;
   
   sep = dot( rsun , rsc )/( rsunmag*rscmag ) ;
   if  sep < 0 
       light = 1 ; % not in any shadow
       sepang = acos( sep ) ; 
       sath = rscmag*cos( sepang ) ;
       satv = rscmag*sin( sepang ) ;
       penv = tan( alphap )*( ( radEarth / sin( alphap ) ) + sath ) ;
       if satv <= penv
           light = 0 ; % penumbra
           umbv = tan( alphau )*( ( radEarth / sin( alphau ) ) - sath ) ;
           if satv <= umbv
               light = 0 ; % umbra
           end
       end
   else 
       light = 1 ; % light and not on a possible shadow side
   end
   
end


function ap = SRPacceleration( r , v , A , m , JDo , t , solarsail) 
% if solar sail is one then it increases the area by 100 when sun and
% velocity are aligned
    p = 4.57e-6 ;
    [ light , s ] = ValladoShadow( r , JDo , t ) ;
    if solarsail == 1 
        alignment = ( dot( v , s )/( norm(s)*norm(v) ) ) ;
        if alignment < -.9
            A = A*100 ;
            t
        end
    end 
    cr = 1.2 ;
    ap = (- ( p * cr * A * s * light )/m)*1e-3 ;
    
end




%% N Body Code
function p = NBod(r,NBodInp,t,jd)
    MassSun = 1.989E30 ; %kg
    GravConst = 6.6743015E-11/1000^3 ; %km^3/(Kg-s^2)
    tdays = t/(24*60*60) ;
    p = 0 ;
    if NBodInp(1) == 1
        [~ , ~ , r_S] = solar_position(jd+tdays) ;
        r_S_SC = (r_S - r) ; % Invert(Sun to Earth - Earth to S/C = Sun to S/C
        p = GravConst*MassSun*(r_S_SC/norm(r_S_SC)^3-r_S/norm(r_S)^3) ;
    end
    if NBodInp(2) == 1
        p = p + 0;
    end
end
%% Drag Code
function p = Dragger(State,Cd,Area,Mass)
v = State(4:6) ;
r = State(1:3) ;
EarthOmeg = [0;0;72.9211E-6] ; % Rad/s
VelAtmo = cross(EarthOmeg,r) ;
VelDiff = (v-VelAtmo) ;
Dens = atmosphere(norm(r)-6378) ; %Finds Density of air
p = -.5*Cd*Area/Mass*Dens*norm(VelDiff*1000)^2*v/norm(v)/1000 ;
end
%% TLE output to RV
function State = TLE2RV(DateVec,Epochtime,mu,MeanMot,MeanAnom,Ecc,Inc,RAAN,ArgOPerigee)
%For initial at EPOCH time
EpochStr = num2str(Epochtime) ;
Epoch.y = str2double(EpochStr(1:2)) ;
Epoch.d = str2double(EpochStr(3:end)) ;
JDate = JDconv(Epoch) ; %in days since beginning

%Chose Date as 10/10/2019 12:00
JDateToday = juliandate(DateVec) ;
TimeSince = (JDateToday - JDate)*24*60*60 ; % Time in seconds since TLE

% M = Mean Anomaly
% N = Mean Motion

SMA = (mu/(MeanMot*2*pi/(24*60*60))^2)^(1/3);
h = sqrt(SMA*(1-Ecc^2)*mu);

% Need TA
if MeanAnom < 2*pi
    E = MeanAnom+Ecc/2 ; %first guess if MA is less than 2*pi
else
    E = MeanAnom-Ecc/2 ; %first guess if MA is more than 2*pi
end
%Use newtons method to find E
E1 = Inf ; %ensures step into while loop
E2 = 1 ; %ensures step into while
while abs(E2 - E1) > 1E-8 %Sets tolerance for newton iteration
    E1 = E+(MeanAnom-E+Ecc*sin(E))/(1-Ecc*cos(E)) ; %first step in newton iteration (looped)
    E2 = E ; %stores E value for the while loop tolerance
    E = E1 ; %Sets the new value of E
end
TAF = 2*atan2(tan(E/2),(sqrt((1-Ecc)/(1+Ecc)))) ; %solves for the final True anomaly
TA = TAF/pi*180 ;%now in degrees


% Initial State
State = Orb2ECI(h,mu,Ecc,TA,Inc,RAAN,ArgOPerigee) ;
 function [JD] = JDconv(epoch)
    JD = 2451545.0+365.25*epoch.y+epoch.d;
 end
function [State] = Orb2ECI(h,mu,ecc,TA,inc,RAAN,ArgOfPerigee)
    % Combines elements to perifocal and perifocal to ECI
    PeriState = PERISTATE(h,mu,ecc,TA) ;
    [State] = PeriToECI(PeriState,inc,RAAN,ArgOfPerigee) ;
end
    %% Perifocal to State
    function State = PERISTATE(h,mu,ecc,TA)
    %input 4 values to get the state vector in a perifocal frame
    State(1:3,1) = h^2/(mu*(1+ecc*cosd(TA)))*[cosd(TA);sind(TA);0] ; %position in perifocal

    State(4:6,1) = mu/h*[-sind(TA);(ecc+cosd(TA));0] ; %Velocity in Perifocal
    end

    %% Perifocal to ECI
    function [State] = PeriToECI(State,inc,RAAN,ArgOfPerigee)
    %state must be input as a vector
        function [Mat] = R3(Angle) % generates a 3rd axis rotation matrix
            Mat = [cosd(Angle),-sind(Angle),0;sind(Angle),cosd(Angle),0;0,0,1] ;
        end
        function [Mat] = R1(Angle) % generates a 1st axis rotation matrix
            Mat = [1,0,0;0,cosd(Angle),-sind(Angle);0,sind(Angle),cosd(Angle)];
        end
    State(1:3) = R3(RAAN)*R1(inc)*R3(ArgOfPerigee)*State(1:3) ; %transforms the position into perifocal
    State(4:6) = R3(RAAN)*R1(inc)*R3(ArgOfPerigee)*State(4:6) ; %transforms the velocity into perifocal
    end
end
%% Curtis Exponential Model For atmosphere
%wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
function density = atmosphere(z)
% ATMOSPHERE calculates density for altitudes from sea level
% through 1000 km using exponential interpolation.
%wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
%...Geometric altitudes (km):

ha = ...
[ 0 25 30 40 50 60 70 ...
80 90 100 110 120 130 140 ...
150 180 200 250 300 350 400 ...
450 500 600 700 800 900 1000];
%...Corresponding densities (kg/m^3) from USSA76:
r = ...
[1.225 4.008e-2 1.841e-2 3.996e-3 1.027e-3 3.097e-4 8.283e-5 ...
1.846e-5 3.416e-6 5.606e-7 9.708e-8 2.222e-8 8.152e-9 3.831e-9 ...
2.076e-9 5.194e-10 2.541e-10 6.073e-11 1.916e-11 7.014e-12 2.803e-12 ...
1.184e-12 5.215e-13 1.137e-13 3.070e-14 1.136e-14 5.759e-15 3.561e-15];
%...Scale heights (km):
H = ...
[ 7.310 6.427 6.546 7.360 8.342 7.583 6.661 ...
5.927 5.533 5.703 6.782 9.973 13.243 16.322 ...
21.652 27.974 34.934 43.342 49.755 54.513 58.019 ...
60.980 65.654 76.377 100.587 147.203 208.020];
%...Handle altitudes outside of the range:
if real(z) > 1000
    z = 1000;
elseif real(z) < 0
    z = 0;
end
%...Determine the interpolation interval:
for j = 1:27
    if z >= ha(j) && z < ha(j+1)
        i = j;
    end
end
if z == 1000
    i = 27;
end

%...Exponential interpolation:
density = r(i)*exp(-(z - ha(i))/H(i)); %atmopshere
end
%% Orb Elements From RV (FOR SMA)
function [ha,inca,ecca,RAANa,ARGOPER,TA,SMA] = orbitvars(pos1,vel1,mu)
% Input position and velocity, and the large body gravity constant (use km/s and km and km^3/s^2)
[M,N] = size(pos1) ;
for oo = 1:M+N-3
    pos = pos1(:,oo) ;
    vel = vel1(:,oo) ;
    dist = norm(pos) ; %gives the distance to the center of the earth
    speed = norm(vel) ; %gives us the speed at given distance 
    vr = dot(pos,vel)/dist ; %gives us the radial velocity
    hvec = cross(pos,vel) ; %gives us the angular momentum vector
    ha(oo) = norm(hvec) ; %gives the magnitude of angular momentum 
    inca(oo) = acosd(hvec(3)/ha(oo)) ; %gives the inclination in degrees
    nodevec = cross([0;0;1],hvec) ; %gives the node line
    node = norm(nodevec) ; %gives mag of node line
    if nodevec(2) > 0
        RAANa(oo) = acosd(nodevec(1)/node) ; %gives us the right ascencion of the ascending node
    else
        RAANa(oo) = 360 - acosd(nodevec(1)/node) ; %gives us the right ascencion of the ascending node
    end
    eccvec = 1/mu*((speed^2-mu/dist)*pos-dist*vr*vel) ; %gives the ecc vector
    ecca(oo) = norm(eccvec) ; %gives the actual ecc value
    if eccvec(3) > 0 
        ARGOPER(oo) = acosd(dot(nodevec,eccvec)/(node*ecca(oo))) ; %gives the arg of perigee in degrees
    else 
        ARGOPER(oo) = 360 - acosd(dot(nodevec,eccvec)/(node*ecca(oo))) ; %gives the arg of perigee in degrees
    end
    if vr >= 0
        TA(oo) = acosd(dot(eccvec,pos)/(ecca(oo)*dist)) ; %gives the true anomaly
    else 
        TA(oo) = 360 - acosd(dot(eccvec,pos)/(ecca(oo)*dist)) ; %gives the true anomaly
    end
    SMA(oo) = ha(oo)^2/(mu*(1-ecca(oo)^2)) ; %gives the semimajor axis
end
end
function [dV1 , dV2 , ecc] = MockHohmannCalculator(R1,V1,R2,V2,mu)
%Input the final and initial radius of two circular orbits, output deltaV's
%of each and ecc of transfer orbit (need mu input)
V10 = V1 ;
V1f = V2 ;
ecc = (R2-R1)/(R2+R1) ; %Ecc of the transfer orbit
h = sqrt(R1*(1+ecc)*mu) ; %Specific angular momentum of transfer orbit (using perigee values)
Vperi = mu/h*(1+ecc) ; %Velocity at transfer Perigee
Vapp = mu/h*(1-ecc) ; %Velocity at transfer Perigee
%Separate burn 1 and burn 2
dV1 = abs(Vperi-V10) ; %km/s dV for burn 1
dV2 = abs(Vapp-V1f) ; %km/s dV for burn 2 (at apogee)
end
%% Event Function
function [value, isterminal, direction] = MyEvent(~, y,~,~,~,~,~,~,~,~,~)
value      =    norm(y(1:3)) - 120 - 6378 ;
isterminal =    1;   % Stop the integration
direction  =    -1;
end
%% Curtis Code For Solar Location
% D.45 Algorithm 12.2: Calculate the geocentric position of the sun at a
% given epoch
% Function file: solar_position.m
% wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
function [lamda , eps , r_S] = solar_position(jd)
%
% This function calculates the geocentric equatorial position vector
% of the sun, given the Julian date.
%
% User M-functions required: None
% -------------------------------------------------------------------------
%...Astronomical unit (km):
AU = 149597870.691;
%...Julian days since J2000:
n = jd - 2451545;
%...Julian centuries since J2000:
cy = n/36525;
%...Mean anomaly (deg{:
M = 357.528 + 0.9856003*n;
M = mod(M,360);
%...Mean longitude (deg):
L = 280.460 + 0.98564736*n;
L = mod(L,360);
%...Apparent ecliptic longitude (deg):
lamda = L + 1.915*sind(M) + 0.020*sind(2*M);
lamda = mod(lamda,360);
%...Obliquity of the ecliptic (deg):
eps = 23.439 - 0.0000004*n;
%...Unit vector from earth to sun:
u = [cosd(lamda); sind(lamda)*cosd(eps); sind(lamda)*sind(eps)];
%...Distance from earth to sun (km):
rS = (1.00014 - 0.01671*cosd(M) - 0.000140*cosd(2*M))*AU;
%...Geocentric position vector (km):
r_S = rS*u;
end %solar_position
%% Lambert's Problem Calculator
function [V1 , V2] = LambertoCalc(R1,R2,TimePass,mu,Direction)
%Uses the Stumpff equations and Lambert's Theorem to take two position
%vectors and time passage, to give velocity vectors at both positions.
%input 1 or -1 for direction
%Can Input muearth, or any other mu
R1n = norm(R1) ; %magnitude of R1
R2n = norm(R2) ; %magnitude of R2
%Finding true anom diff
Crossed = cross(R1,R2) ; %Used to check short and long
if Crossed(3) >= 0
dThetaPro = acosd(dot(R1,R2)/(norm(R1)*norm(R2))) ; %difference in
% true anomaly (for prograde)
dThetaRetro = 360 - acosd(dot(R1,R2)/(norm(R1)*norm(R2))) ; %for retrograde
else
dThetaPro = 360 - acosd(dot(R1,R2)/(norm(R1)*norm(R2))) ; %difference in true anomaly
dThetaRetro = acosd(dot(R1,R2)/(norm(R1)*norm(R2))) ; %Retrograde
end
%Selection of Prograde or Retrograde
if Direction == -1
dTheta = dThetaRetro ;
else
dTheta = dThetaPro ;
end
%Functions for Lamberts
VFe = @(g,gdot,R0,RF) 1/g*(gdot*RF-R0) ; % For the final velocity,
% where R0 and V0 are vectors
V0e = @(g,f,R0,RF) 1/g*(RF-f*R0) ; %for initial velocity vector
ge = @(A,y) A*sqrt(y/mu) ; %for g (when A and y is known)
fe = @(y) 1 - y/R1n ; %for f when y and R0 are known
gdote = @(y)(1-y/R2n) ; %for gdot with y and RF known
%A and y and Chi
%Use Stumpf equations (see last functions
Ae = @(dTheta) sqrt(R1n*R2n)*sind(dTheta)/(sqrt(1-cosd(dTheta))) ; %For 'A'
ye = @(A,Z) R1n + R2n + A*(Z*StumS(Z)-1)/sqrt(StumC(Z)) ; %for y
Chi = @(y,Z) sqrt(y)/sqrt(StumC(Z)) ; %Finding Chi includes both y and z
%Iterate to find Z (Bisection Method)
%Guesses'
oo = 0 ;
Z = 0 ; %initial guesses
Zupper = 4*pi^2 ; Zlower = -4*pi^2 ; %Z bounds for bisection
guessT = Inf ; %ensures stepping into while loop
while abs(real(guessT) - TimePass) > 10E-8 %To find Z where time equals 20 minutes
    oo = oo + 1;
%Caluclate A, y, and Chi, most use stumpf eq's
AN = Ae(dTheta) ; %for A
YN= ye(AN,Z) ; %Z will change and for Y
ChiN = Chi(YN,Z) ; %for Chi
%For new T calculated
guessT = ChiN^3*StumS(Z)/sqrt(mu)+AN*sqrt(YN)/sqrt(mu) ; %fordelta T
if guessT < TimePass %bisection
Zlower = Z ; %Shrinks bounds
end
if guessT > TimePass %bisection
Zupper = Z ; %shrinks bounds
end
if guessT ==TimePass
Zupper = Z ; %Makes the bounds the same
Zlower = Z ; %Because Z is exact
end
Z = (Zupper+Zlower)/2 ; %new Z value
if oo > 1000
    break
end
end
%Use Functions Made Above to Find V1 and V2
%Recalculate A, y, and Chi
%For Prograde (slot 1 in Z)
A1 = Ae(dTheta) ; %for A
Y1 = ye(A1,Z) ; %using iterated Z to find Y
%Finding F,Fdot,G,Gdot
g1 = ge(A1,Y1) ; %for g
gdot1 = gdote(Y1) ; %for gdot
f1 = fe(Y1) ; %for f
%Finding Velocities
if oo > 10000
    V1 = [Inf;Inf;Inf] ; V2 = [Inf;Inf;Inf] ;
else
V1 = V0e(g1,f1,R1,R2) ; %Finding V1 prograde velocity vector
V2 = VFe(g1,gdot1,R1,R2) ;%Finding V2 prograde velocity vector
end
%Stumpf Equations Used Above
%Stumpf Equations-Given By lecture
%Separated Functions, Separated Outputs
function [C] = StumC(Z)
if Z > 0
C = (1-cos(sqrt(Z)))/Z ;
end
if Z < 0
C = (cosh(sqrt(-Z))-1)/(-Z) ;
end
if Z == 0
C = 1/2 ;
end
end
function [S] = StumS(Z)
if Z > 0
S = (sqrt(Z) - sin(sqrt(Z)))/Z^1.5 ; %Stumph Eq's
end
if Z < 0
S = (sinh(sqrt(-Z))-sqrt(-Z))/(-Z)^1.5 ; %if Z is negative
end
if Z == 0
S = 1/6 ;
end
end
end
function GlobeTrotter()
grs80 = referenceEllipsoid('grs80','km');
domeRadius =  6378;  % km
domeLat =  39;       % degrees
domeLon = -77;       % degrees
domeAlt = 0;         % km
[x,y,z] = sphere(20);
xEast  = domeRadius * x;
yNorth = domeRadius * y;
zUp    = domeRadius * z;
zUp(zUp < 0) = 0;
axis equal
view(60,60)
axis equal
ax = axesm('globe','Geoid',grs80,'Grid','off', ...
    'GLineWidth',1,'GLineStyle','-',...
    'Gcolor',[0.9 0.9 0.1],'Galtitude',100);
ax.Position = [0 0 1 1];
load topo
geoshow(topo,topolegend,'DisplayType','texturemap','HandleVisibility','off')
demcmap(topo)
land = shaperead('landareas','UseGeoCoords',true);
plotm([land.Lat],[land.Lon],'Color','black','HandleVisibility','off')
rivers = shaperead('worldrivers','UseGeoCoords',true);
plotm([rivers.Lat],[rivers.Lon],'Color','blue','HandleVisibility','off')
hold on
end



