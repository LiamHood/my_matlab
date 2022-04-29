close all; clear; clc;

T = 5;
V = 15;
dt = .1;
tend = 100;

t_tes(1) = 0;
x_tes(1) = -80;
y_tes(1) = 0;
theta_tes(1) = 0;
thetadot1(1) = 0;

ii = 1;
while t_tes(end) < 25
    [t_tes(ii+1), x_tes(ii+1), y_tes(ii+1), theta_tes(ii+1)] = prop_desired(t_tes(ii), x_tes(ii), y_tes(ii), theta_tes(ii), dt, V);
    ii = ii + 1;
end

t1(1) = 0;
x1(1) = -80;
y1(1) = 0;
theta1(1) = 0;
thetadot1(1) = 0;
t2(1) = 0;
x2(1) = -80;
y2(1) = 0;
theta2(1) = pi;
thetadot2(1) = 0;
ii = 1;
while t1(end) < tend
    r1(ii) = .4 + (pi/2)*sin((2*pi/T)*t1(ii));
    [t1(ii+1), x1(ii+1), y1(ii+1), theta1(ii+1), thetadot1(ii+1)] = prop(t1(ii), x1(ii), y1(ii), theta1(ii), thetadot1(ii), r1(ii), dt, T, V);
    r2(ii) = -(.3 + (pi/2)*cos((2*pi/T)*t2(ii)));
    [t2(ii+1), x2(ii+1), y2(ii+1), theta2(ii+1), thetadot2(ii+1)] = prop(t2(ii), x2(ii), y2(ii), theta2(ii), thetadot2(ii), r2(ii), dt, T, V);
    ii = ii + 1;
end

tra_input =    [theta1(2:end), theta2(2:end);
                theta1(1:end-1), theta2(1:end-1);
                thetadot1(1:end-1), thetadot2(1:end-1);
                dt*ones(1,length(r1)), dt*ones(1,length(r2))];
tra_target = [r1, r2];

% net = feedforwardnet(2,'trainlm');
% net = configure(net, tra_input, tra_target);
% net.layers{1}.transferFcn = 'tansig';
% % net.layers{2}.transferFcn = 'tansig';
% net.trainParam.show = 50;
% net.trainParam.lr = .05;
% net.trainParam.epochs = 10000;
% net.trainParam.goal = 1e-12;
% net.divideParam.trainRatio = .750;
% net.divideParam.valRatio = .250;
% net.divideParam.testRatio = .00;
% net = train(net, tra_input, tra_target);
% disp(perform(net,tra_input, tra_target))

load("net_tenth_second.mat")
% load("net_1_second.mat")

net_tra_output = net(tra_input);
nr1 = net_tra_output(1:end/2);
nr2 = net_tra_output(end/2+1:end);
nt1(1) = 0;
nx1(1) = -80;
ny1(1) = 0;
ntheta1(1) = 0;
nthetadot1(1) = 0;
nt2(1) = 0;
nx2(1) = -80;
ny2(1) = 0;
ntheta2(1) = pi;
nthetadot2(1) = 0;
ii = 1;
while nt1(end) < tend
    [nt1(ii+1), nx1(ii+1), ny1(ii+1), ntheta1(ii+1), nthetadot1(ii+1)] = prop(nt1(ii), nx1(ii), ny1(ii), ntheta1(ii), nthetadot1(ii), nr1(ii), dt, T, V);
    [nt2(ii+1), nx2(ii+1), ny2(ii+1), ntheta2(ii+1), nthetadot2(ii+1)] = prop(nt2(ii), nx2(ii), ny2(ii), ntheta2(ii), nthetadot2(ii), nr2(ii), dt, T, V);
    ii = ii + 1;
end

ntt1(1) = 0;
nxt1(1) = -80;
nyt1(1) = 0;
nthetat1(1) = 0;
nthetadott1(1) = 0;
ntt2(1) = 0;
nxt2(1) = -80;
nyt2(1) = 0;
nthetat2(1) = 0;
nthetadott2(1) = 0;
ntt3(1) = 0;
nxt3(1) = -80;
nyt3(1) = 0;
nthetat3(1) = 0;
nthetadott3(1) = 0;
ii = 1;
while ntt1(end) < 25
    % Calculated based on needed heading angles, propagating normally 
    theta_needed(ii) = atan((x_tes(ii+1)-nx1(ii))/(y_tes(ii+1)-ny1(ii)));
    if ii > 1
        dif = abs(theta_needed(ii-1) - theta_needed(ii));
        if dif > pi
            theta_needed(ii) = theta_needed(ii) + 2*pi;
        end
    end
    theta_calc = atan((x_tes(ii+1)-x_tes(ii))/(y_tes(ii+1)-y_tes(ii)));
    nrt1(ii) = net([theta_needed(ii); nthetat1(ii); nthetadott1(ii); dt]);
%     disp(ii)
%     disp("to nominal")
%     disp(theta_needed(ii))
%     disp(theta_calc)
    [ntt1(ii+1), nxt1(ii+1), nyt1(ii+1), nthetat1(ii+1), nthetadott1(ii+1)] = prop(ntt1(ii), nxt1(ii), nyt1(ii), nthetat1(ii), nthetadott1(ii), nrt1(ii), dt, T, V);
    
    % Calculated based on current heading angle and original desired
    nrt2(ii) = net([theta_tes(ii+1); nthetat2(ii); nthetadott2(ii); dt]);
%     disp("heading")
%     disp(theta_tes(ii+1))
%     disp(nrt2(ii))
    [ntt2(ii+1), nxt2(ii+1), nyt2(ii+1), nthetat2(ii+1), nthetadott2(ii+1)] = prop(ntt2(ii), nxt2(ii), nyt2(ii), nthetat2(ii), nthetadott2(ii), nrt2(ii), dt, T, V);
    
    % One step ahead controller
    nrt3(ii) = net([theta_tes(ii+1); theta_tes(ii); nthetadott3(ii); dt]);
    [ntt3(ii+1), nxt3(ii+1), nyt3(ii+1), nthetat3(ii+1), nthetadott3(ii+1)] = prop(t_tes(ii), x_tes(ii), y_tes(ii), theta_tes(ii), nthetadott3(ii), nrt3(ii), dt, T, V);
    ii = ii + 1;
end

figure(1)
axis([-100, 100, -100, 100],'square')
hold on
plot(x_tes,y_tes,"k")
plot(x1, y1, "g")
plot(x2, y2, "b")
plot(nx1, ny1, ".r")
plot(nx2, ny2, ".r")
plot(nxt1, nyt1, "-.")
plot(nxt2, nyt2, "*")
plot(nxt3, nyt3, ".")
xlabel("x [m]")
ylabel("y [m]")
title("Path")
legend("Target", "Training1", "Training2", "Net Training 1", "Net Training 2", ...
    "Propagation with Control based on Nominals" , "Propagation with Control based on Actuals", "One Step Ahead", "Location", "bestoutside")
%     "Training1", "Training2", "Net Training 1", "Net Training 2", ...
    

figure(2)
hold on
plot(t_tes, theta_tes, "k")
plot(t1, theta1, "g")
plot(t2, theta2, "b")
plot(nt1, ntheta1, ".r")
plot(nt2, ntheta2, ".r")
plot(ntt1, nthetat1, ".")
plot(ntt2, nthetat2, "*")
plot(ntt3, nthetat3, ".")
xlabel("Time [s]")
ylabel("Heading Angle [radians]")
title("Heading Angle over Time")
legend("Target", "Training1", "Training2", "Net Training 1", "Net Training 2", ...
    "Propagation with Control based on Nominals" , "Propagation with Control based on Actuals", "One Step Ahead", "Location", "bestoutside")

figure(3)
hold on
plot(t1, thetadot1, "g")
plot(t2, thetadot2, "b")
plot(nt1, nthetadot1, ".r")
plot(nt2, nthetadot2, ".r")
plot(ntt1, nthetadott1, ".")
plot(ntt2, nthetadott2, "*")
plot(ntt3, nthetadott3, ".")
xlabel("Time [s]")
ylabel("Turning Rate [radians/s]")
title("Turning Rate over Time")
legend("Training1", "Training2", "Net Training 1", "Net Training 2", ...
    "Propagation with Control based on Nominals" , "Propagation with Control based on Actuals", "One Step Ahead", "Location", "bestoutside")

figure(4)
hold on
plot(t1(1:end-1), r1, "g")
plot(t2(1:end-1), r2, "b")
plot(nt1(1:end-1), nr1, ".r")
plot(nt2(1:end-1), nr2, ".r")
plot(ntt1(1:end-1), nrt1, ".")
plot(ntt2(1:end-1), nrt2, "*")
plot(ntt3(1:end-1), nrt3, ".")
xlabel("Time [s]")
ylabel("Rudder Angle [radians]")
title("Rudder Command over Time")
legend("Training1", "Training2", "Net Training 1", "Net Training 2", ...
    "Propagation with Control based on Nominals" , "Propagation with Control based on Actuals", "One Step Ahead", "Location", "bestoutside")


function [t_next, x_next, y_next, theta_next, thetadot_next] = prop(t, x, y, theta, thetadot, r, dt, T, V)
    thetadot_next = thetadot + dt*((r - thetadot)/T);
    theta_next = theta + dt*((thetadot_next + thetadot)/2);
    x_next = x + dt*V*sin((theta_next + theta)/2);
    y_next = y + dt*V*cos((theta_next + theta)/2);
    t_next = t + dt;
end

function [t_next, x_next, y_next, theta_next] = prop_desired(t, x, y, theta, dt, V)
    theta_next = theta + dt*(.2+.2*theta/(2*pi));
    x_next = x + dt*V*sin((theta_next + theta)/2);
    y_next = y + dt*V*cos((theta_next + theta)/2);
    t_next = t + dt;
end