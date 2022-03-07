clear; close all; clc;

%% Constants
h = 1;
l1 = 1.6;
l2 = 1.6;
n = 5;
m = 5;

x2math = @(theta1, theta2, l1, l2) l1*cosd(theta1 - 90) + l2*sind(theta2 + theta1 - 180);
y2math = @(theta1, theta2, l1, l2, h) h + l1*sind(theta1 - 90) - l2*cosd(theta2 + theta1 - 180);

%% Training Data
min1 = 105;
max1 = 165;
min2 = 35;
max2 = 90;
t_th1 = linspace(min1,max1,n)';
t_th2 = linspace(min2,max2,m)';
kk = 0;
for ii = 1:n
    for jj = 1:m
        if (t_th1(ii) + t_th2(jj)) < 220 && (t_th1(ii) + t_th2(jj)) > 155
            kk = kk + 1;
%             n*(ii-1)+jj
            t_x2(kk,1) = x2math(t_th1(ii), t_th2(jj), l1, l2);
            t_y2(kk,1) = y2math(t_th1(ii), t_th2(jj), l1, l2, h);
            tlra_tra(kk,:) = [t_x2(kk,1), t_y2(kk,1), t_th1(ii), t_th2(jj)];
        end
    end
end
figure
axis([-1,3,-1,3],'square')
hold on
plot(t_x2, t_y2, '.')
plot([1,2], [0,0], "r")
plot([1,2], [1,1], "r")
plot([1,1], [0,1], "r")
plot([2,2], [0,1], "r")
plot([0,0], [0,1], "k")
hold off
title("Training Data")
xlabel("x1")
ylabel("x2")


%% Testing Data
test_x1 = linspace(1,2,11)';
test_y1_1 = test_x1 - 1;
test_y1_2 = 2 - test_x1;
theta = (0:5:355)';
tl = length(theta);
test_x2 = 1.5 + 0.3*cosd(theta);
test_y2 = 0.5 + 0.3*sind(theta);
test_x3 = 1.5 + 0.6*cosd(theta);
test_y3 = 0.5 + 0.6*sind(theta);

tlra_tes = [test_x1, test_y1_1, zeros(11,1), zeros(11,1);
            test_x1, test_y1_2, zeros(11,1), zeros(11,1);
            test_x2, test_y2, zeros(tl,1), zeros(tl,1);
            test_x3, test_y3, zeros(tl,1), zeros(tl,1)];

% figure
% % axis([1,2,0,1],'square')
% axis([-1,3,-1,3],'square')
% hold on
% plot(t_x2, t_y2, '.')
% plot(test_x1, test_y1_1, "o")
% plot(test_x1, test_y1_2, "o")
% plot([1,2], [0,0], "r")
% plot([1,2], [1,1], "r")
% plot([1,1], [0,1], "r")
% plot([2,2], [0,1], "r")
% plot([0,0], [0,1], "k")
% hold off
% title("Test 1")
% xlabel("x1")
% ylabel("x2")
% 
% figure
% axis([-1,3,-1,3],'square')
% hold on
% plot(t_x2, t_y2, '.')
% plot(test_x2, test_y2, "o")
% plot([1,2], [0,0], "r")
% plot([1,2], [1,1], "r")
% plot([1,1], [0,1], "r")
% plot([2,2], [0,1], "r")
% plot([0,0], [0,1], "k")
% hold off
% title("Test 2")
% xlabel("x1")
% ylabel("x2")
% 
% figure
% axis([-1,3,-1,3],'square')
% hold on
% plot(t_x2, t_y2, '.')
% plot(test_x3, test_y3, "o")
% plot([1,2], [0,0], "r")
% plot([1,2], [1,1], "r")
% plot([1,1], [0,1], "r")
% plot([2,2], [0,1], "r")
% plot([0,0], [0,1], "k")
% hold off
% title("Test 3")
% xlabel("x1")
% ylabel("x2")

%% Show Test

load("test_9.mat")

r_x2(:,1) = x2math(results(:,1), results(:,2), l1, l2);
r_y2(:,1) = y2math(results(:,1), results(:,2), l1, l2, h);
r1_x = r_x2(1:22);
r1_y = r_y2(1:22);
r2_x = r_x2(23:94);
r2_y = r_y2(23:94);
r3_x = r_x2(95:end);
r3_y = r_y2(95:end);
x1err = r1_x - [test_x1; test_x1];
y1err = r1_y - [test_y1_1; test_y1_2];
rms1 = sqrt(mean([x1err;y1err].^2));
fprintf("RMS for X test is :%f \n", rms1)

x2err = r2_x - test_x2;
y2err = r2_y - test_y2;
rms2 = sqrt(mean([x2err;y2err].^2));
fprintf("RMS for small circle test is :%f \n", rms2)

x3err = r3_x - test_x3;
y3err = r3_y - test_y3;
rms3 = sqrt(mean([x3err;y3err].^2));
fprintf("RMS for large circle test is :%f \n", rms3)


rms_all = sqrt(mean([x1err;y1err;x2err;y2err;x3err;y3err].^2));
fprintf("RMS for all tests is :%f \n", rms_all)

figure
axis([0,2.5,-.5,2],'square')
hold on
plot(test_x1, test_y1_1, "o")
plot(test_x1, test_y1_2, "o")
plot(r1_x, r1_y, "*")
plot([1,2], [0,0], "r")
plot([1,2], [1,1], "r")
plot([1,1], [0,1], "r")
plot([2,2], [0,1], "r")
plot([0,0], [0,1], "k")
hold off
title("Results of Test 1")
xlabel("x1")
ylabel("x2")
legend("Target Points", "Robot position based on NN output")

figure
axis([0,2.5,-.5,2],'square')
hold on
plot(test_x2, test_y2, "o")
plot(r2_x, r2_y, "*")
plot([1,2], [0,0], "r")
plot([1,2], [1,1], "r")
plot([1,1], [0,1], "r")
plot([2,2], [0,1], "r")
plot([0,0], [0,1], "k")
hold off
title("Results of Test 2")
xlabel("x1")
ylabel("x2")
legend("Target Points", "Robot position based on NN output")


figure
axis([0,2.5,-.5,2],'square')
hold on
plot(test_x3, test_y3, "o")
plot(r3_x, r3_y, "*")
plot([1,2], [0,0], "r")
plot([1,2], [1,1], "r")
plot([1,1], [0,1], "r")
plot([2,2], [0,1], "r")
plot([0,0], [0,1], "k")
hold off
title("Results of Test 3")
xlabel("x1")
ylabel("x2")
legend("Target Points", "Robot position based on NN output")
