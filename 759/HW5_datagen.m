clear; close all; clc;

theta = @(t,J,T) (T/(2*J))*t^2;
omega = @(t,J,T) (T/J)*t;
J = 1.5;
l = 2;

t = 0:.1:1.1;

Ttra = [1,2,3,4];

Ttes = [.5, 1.5, 2.5, 3.5, 4.5, 5.5];

n = length(Ttra);
m = length(t);
for ii = 1:n
    for jj = 1:m
        tra((ii-1)*m+jj,:) = [t(jj),theta(t(jj),J,Ttra(ii)),omega(t(jj),J, ...
            Ttra(ii)),Ttra(ii),l*cos(theta(t(jj),J,Ttra(ii))),l*sin(theta(t(jj),J,Ttra(ii)))];
    end
end

n = length(Ttes);
m = length(t);
for ii = 1:n
    for jj = 1:m
        tes((ii-1)*m+jj,:) = [t(jj),theta(t(jj),J,Ttes(ii)),omega(t(jj),J, ...
            Ttes(ii)),Ttes(ii),l*cos(theta(t(jj),J,Ttes(ii))),l*sin(theta(t(jj),J,Ttes(ii)))];
    end
end

figure
axis([-3,3,-3,3],'square')
hold on
for ii = 1:4
    plot(tra(1+n*(ii-1):n*ii,5),tra(1+n*(ii-1):n*ii,6),'.')
end
xlabel('x1')
ylabel('x2')
title('Training Data')
legend('T=1','T=2','T=3','T=4','Location','southwest')

figure
axis([-3,3,-3,3],'square')
hold on
for ii = 1:m
    plot(tes(1+n*(ii-1):n*ii,5),tes(1+n*(ii-1):n*ii,6),'.')
end
xlabel('x1')
ylabel('x2')
title('Testing Data')
legend('T=.5','T=1.5','T=2.5','T=3.5','T=4.5','T=5.5','Location','southwest')

clear; close all; clc;
load("t_4_ESin.mat")

theta = @(t,J,T) (T/(2*J))*t^2;
omega = @(t,J,T) (T/J)*t;
J = 1.5;
l = 2;

t = 0:.1:1.1;

Ttra = [1,2,3,4];

Ttes = [.5, 1.5, 2.5, 3.5, 4.5, 5.5];

n = length(Ttra);
m = length(t)-1;
for ii = 1:n
    for jj = 1:m
        tra((ii-1)*m+jj,:) = [Ttra(ii),...
            theta(t(jj),J,Ttra(ii)),...
            omega(t(jj),J,Ttra(ii)),...
            l*cos(theta(t(jj+1),J,Ttra(ii))),...
            l*sin(theta(t(jj+1),J,Ttra(ii))),...
            theta(t(jj+1),J,Ttra(ii)),...
            omega(t(jj+1),J,Ttra(ii))];
    end
end

n = length(Ttes);
for ii = 1:n
    for jj = 1:m
        tes((ii-1)*m+jj,:) = [Ttes(ii),...
            theta(t(jj),J,Ttes(ii)),...
            omega(t(jj),J,Ttes(ii)),...
            l*cos(theta(t(jj+1),J,Ttes(ii))),...
            l*sin(theta(t(jj+1),J,Ttes(ii))),...
            theta(t(jj+1),J,Ttes(ii)),...
            omega(t(jj+1),J,Ttes(ii))];
    end
end

% n = length(Ttra);
% m = length(t);
% for ii = 1:n
%     for jj = 1:m
%         tra((ii-1)*m+jj,:) = [Ttra(ii),t(jj),theta(t(jj),J,Ttra(ii)),omega(t(jj),J, ...
%             Ttra(ii)),l*cos(theta(t(jj),J,Ttra(ii))),l*sin(theta(t(jj),J,Ttra(ii)))];
%     end
% end
% 
% n = length(Ttes);
% m = length(t);
% for ii = 1:n
%     for jj = 1:m
%         tes((ii-1)*m+jj,:) = [Ttes(ii), t(jj), theta(t(jj),J,Ttes(ii)), omega(t(jj),J, ...
%             Ttes(ii)), l*cos(theta(t(jj),J,Ttes(ii))), l*sin(theta(t(jj),J,Ttes(ii)))];
%     end
% end


% for ii = 1:n
% figure
% axis([-3,3,-3,3],'square')
% hold on
% plot(results(1+m*(ii-1):m*ii,5),results(1+m*(ii-1):m*ii,6),"*")
% plot(results(1+m*(ii-1):m*ii,1),results(1+m*(ii-1):m*ii,2),".")
% xlabel('x1')
% ylabel('x2')
% tes_title = ['Test at Torque = ', num2str(Ttes(ii))];
% title(tes_title)
% legend('Neural','Math Model','Location','southwest')
% end

load("t_4_ESin.mat")
for ii = 1:n
figure
axis([-3,3,-3,3],'square')
hold on
plot(results(1+m*(ii-1):m*ii,5),results(1+m*(ii-1):m*ii,6),"*r")
plot(results(1+m*(ii-1):m*ii,1),results(1+m*(ii-1):m*ii,2),".b")
xlabel('x1')
ylabel('x2')
tes_title = ['Test at Torque = ', num2str(Ttes(ii))];
title(tes_title)
legend('Neural','Math Model','Location','southwest')
end

for ii = 1:4
figure
axis([-3,3,-3,3],'square')
hold on
plot(tra(1+m*(ii-1):m*ii,4),tra(1+m*(ii-1):m*ii,5),'.')
xlabel('x1')
ylabel('x2')
tes_title = ['Training at Torque = ', num2str(Ttra(ii))];
title(tes_title)
end

% figure
% axis([-3,3,-3,3],'square')
% hold on
% % plot(results(:,1),results(:,2),"*")
% for ii = 1:4
%     plot(tra(1+m*(ii-1):m*ii,4),tra(1+m*(ii-1):m*ii,5),'.')
% end
% hold off
% xlabel('x1')
% ylabel('x2')
% title('Training Data')
% legend('Neural','T=1','T=2','T=3','T=4','Location','southwest')
% 
% figure
% axis([-3,3,-3,3],'square')
% hold on
% % plot(results(:,5),results(:,6),"*")
% for ii = 1:n
%     plot(tes(1+m*(ii-1):m*ii,4),tes(1+m*(ii-1):m*ii,5),'.')
% end
% xlabel('x1')
% ylabel('x2')
% title('Testing Data')
% legend('Neural','T=.5','T=1.5','T=2.5','T=3.5','T=4.5','T=5.5','Location','southwest')
>>>>>>> 04d27428a025975f3a5ea7483208d22b363f9248
