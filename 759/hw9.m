close all; clear; clc;

T = 5;
V = 15;
dt = .1;

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
ii = 1;
while t1(end) < T
    r1(ii) = (pi/2)*sin((2*pi/T)*t1(ii));
    [t1(ii+1), x1(ii+1), y1(ii+1), theta1(ii+1), thetadot1(ii+1)] = prop(t1(ii), x1(ii), y1(ii), theta1(ii), thetadot1(ii), r1(ii), dt, T, V);
    ii = ii + 1;
end

t2(1) = 0;
x2(1) = 0;
y2(1) = 0;
theta2(1) = pi;
thetadot2(1) = pi/4;
ii = 1;
while t2(end) < T
    r2(ii) = -(pi/2)*sin((2*pi/T)*t2(ii));
    [t2(ii+1), x2(ii+1), y2(ii+1), theta2(ii+1), thetadot2(ii+1)] = prop(t2(ii), x2(ii), y2(ii), theta2(ii), thetadot2(ii), r2(ii), dt, T, V);
    ii = ii + 1;
end

tra_input1 = [x1(1:end-1);y1(1:end-1);theta1(1:end-1);thetadot1(1:end-1);dt*ones(1,length(r1))];
tra_output1 = [x1(2:end);y1(2:end);theta1(2:end);thetadot1(2:end);r1];
tra_input2 = [x2(1:end-1);y2(1:end-1);theta2(1:end-1);thetadot2(1:end-1);dt*ones(1,length(r2))];

net1 = feedforwardnet([10,10],'trainlm');
net1 = train(net1, tra_input1, tra_output1);
disp(perform(net1,tra_input1, tra_output1))
net1_output1 = net1(tra_input2);
xn1 = net1_output1(1,:);
yn1 = net1_output1(2,:);
thetan1 = net1_output1(3,:);
thetadotn1 = net1_output1(4,:);
rn1 = net1_output1(5,:);
view(net1)

figure(1)
axis([-100, 100, -100, 100],'square')
hold on
plot(x_tes,y_tes,"k")
plot(x1, y1, "g")
plot(x2, y2, "b")
plot(xn1, yn1, "r")
xlabel("x [m]")
ylabel("y [m]")
title("Path")
legend("Target", "Training1", "Training2", "Net 1")

figure(2)
hold on
plot(t_tes, theta_tes, "k")
plot(t1, theta1, "g")
plot(t2, theta2, "b")
plot(t2(2:end), thetan1, "r")
xlabel("Time [s]")
ylabel("Heading Angle [radians]")
title("Heading Angle over Time")
legend("Target", "Training1", "Training2", "Net 1", "Location","southeast")

figure(3)
hold on
plot(t1, thetadot1, "g")
plot(t2, thetadot2, "b")
plot(t2(2:end), thetadotn1, "r")
xlabel("Time [s]")
ylabel("Turning Rate [radians/s]")
title("Turning Rate over Time")
legend("Training1", "Training2", "Net 1", "Location","southeast")

figure(4)
hold on
plot(t1(1:end-1), r1, "g")
plot(t2(1:end-1), r2, "b")
plot(t2(1:end-1), rn1, "r")
xlabel("Time [s]")
ylabel("Rudder Angle [radians]")
title("Rudder Command over Time")
legend("Training1", "Training2", "Net 1", "Location","southeast")


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