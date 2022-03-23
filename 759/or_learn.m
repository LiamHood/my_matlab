clear; clc; close all;
% Set training data
I1 = [0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1];
I2 = [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0];
t3 = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

% Set starting values
w31init = .1;
w32init = .1;
threshold = .5;
c = .25;

w31(1) = w31init;
w32(1) = w32init;
diff = 1;
ii = 1;
while diff > .1 || ii < 18
    x1(ii) = linear_neuron(I1(ii));
    x2(ii) = linear_neuron(I2(ii));
    
    u3(ii) = sum_neuron(x1(ii), w31(ii), x2(ii), w32(ii));
    x3(ii) = step_neuron(u3(ii), threshold);
    err(ii) = t3(ii) - x3(ii);
    [w31(ii+1), w32(ii+1)] = learn(x1(ii), x2(ii), x3(ii), t3(ii), w31(ii), w32(ii), c);
    diff = abs(w31(ii+1)-w31(ii))+abs(w32(ii+1)-w32(ii));
    ii = ii + 1;
end
figure
hold on
plot(1:length(w31), w31)
plot(1:length(w32), w32)
hold off

figure
plot(err)

%% functions

function [w31new, w32new] = learn(x1, x2, x3, t3, w31, w32, c)
    w31new = w31 + c*(t3 - x3)*x1;
    w32new = w32 + c*(t3 - x3)*x2;
end

function x = step_neuron(u, threshold)
    if u > threshold
        x = 1;
    else
        x = 0;
    end
end

function u = sum_neuron(x1, w31, x2, w32)
    u = x1*w31 + x2*w32;
end

function x = linear_neuron(u)
    g = 1;
    x = u*g;
end

