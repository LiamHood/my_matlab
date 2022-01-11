% script written for demo to class

% clear all variables, clear command window, close all objects
clear all; clc; close all

% define Cd0
Cd0 = 0.05;

% define CL and K
CL = [0:0.1:1.5];
K = 0.4;

% calculate CD
for ii=1:length(CL)
    CD(ii) = Cd0 + K*CL(ii)^2;
end

% make plot
figure
plot(CD,CL,'s--')
grid on
ylabel('C_L')
xlabel('C_D')
title(['C_L as f(C_D) for K = ', K])

%% Redo the same calculations without using a for loop
clear
% define Cd0
Cd0 = 0.05;

% define CL and K
CL = [0:0.1:1.5];
K = 0.4;

% calculate CD
CD = Cd0 + K*CL.^2;

% make plot
figure
plot(CD,CL,'rs--')
grid on
ylabel('C_L')
xlabel('C_D')
title(['Same as previous plot, but done without a for loop, K = ', num2str(K)])
% define the angles of attack for which CL was defined
CLalpha = 10;

alpha = CL*CLalpha;

% make plot
figure
plot(alpha, CL./CD,'bs--')
grid on
ylabel('L/D')
xlabel('\alpha')
title(['L/D as f(\alpha), K = ',num2str(K)])

%% clear everything and redo L/D calculations for several K values
% clear

% define Cd0
Cd0 = 0.05;

% define CL and K
jj = 1:length(K);
CL = [0:0.3:1.5];
K =  [0.5:0.1:1];


% calculate CD
CD = Cd0 + K(jj)*CL.^2;

% make plot
figure
hold on
plot(CD,CL','s--')
grid on
ylabel('C_L')
xlabel('C_D')
title(...
    ['Same as previous plot, but done without a for loop and K is now a variable'])

% build legend
for ii = 1:length(K)
    legendEntry{ii} = ['K = ',num2str(K(ii))];
end
legend(legendEntry)

% define the angles of attack for which CL was defined
CLalpha = 10;

alpha = CL*CLalpha;

% make plot
hold on
plot(alpha, CL./CD,'s--')
grid on
ylabel('L/D')
xlabel('\alpha')
title(['L/D as f(\alpha) for various K valus'])


%% learn how to step into a function.  What must you send into this...
% funciton to be done?
myAssignment = isover(2); % ??