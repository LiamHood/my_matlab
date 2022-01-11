clear; close all; clc;
load('flight_test_data.mat')


% filter_747(data_747)
% filter_747(Data7472)
% 
filter_cessna(DataCessna75,"Cessna 75 and 130 knots filtered to values within tolerance")
filter_cessna(DataCessna100,"Cessna 100 knots filtered to values within tolerance")
filter_747(Data747auto,"747 filtered to values within tolerance")


n_ce75 = find(DataCessna75.real_time >= 390 & DataCessna75.real_time <= 390.2, 1, 'first');
n_ce100 = find(DataCessna100.real_time >= 953.8 & DataCessna100.real_time <= 954, 1, 'first');
n_ce130 = find(DataCessna75.real_time >= 1208.4 & DataCessna75.real_time <= 1208.6, 1, 'first');
point(1,1) = plot_point(DataCessna75,n_ce75,73,77,'Cessna 75 knots');
point(2,1) = plot_point(DataCessna100,n_ce100,98,102,'Cessna 100 knots');
point(3,1) = plot_point(DataCessna75,n_ce130,128,132,'Cessna 130 knots');


n_b350 = find(Data747auto.real_time >= 2815.7 & Data747auto.real_time <= 2815.85, 1, 'first');
n_b399 = find(Data747auto.real_time >= 3047.83 & Data747auto.real_time <= 3048, 1, 'first');
n_b450 = find(Data747auto.real_time >= 3178.28 & Data747auto.real_time <= 3178.4, 1, 'first');
point(4,1) = plot_point(Data747auto,n_b350,348,352,'747 350 knots');
point(5,1) = plot_point(Data747auto,n_b399,397,401,'747 399 knots');
point(6,1) = plot_point(Data747auto,n_b450,448,452,'747 450 knots');

results = struct2table(point);
disp(results)
%% Functions

function results = plot_point(data,n,mins,maxs,name)
out.time = data.real_time((n-100):n)-data.real_time(n-100);
out.ba = data.roll__deg((n-100):n);
out.sinkrate = data.VVI__fpm((n-100):n);
out.ta = data.Vtrue_ktas((n-100):n);

figure('Position', [10 10 600 800])
sgtitle(name) 
subplot(3,1,1)
hold on
plot(out.time, out.ba)
plot([0,10],[-5,-5])
plot([0,10],[5,5])
xlabel('Time (s)')
ylabel('Bank Angle (degrees)')
hold off

subplot(3,1,2)
hold on
plot(out.time, out.sinkrate)
plot([0,10],[-500,-500])
plot([0,10],[500,500])
xlabel('Time (s)')
ylabel('Sink Rate (Feet/minute)')
hold off

subplot(3,1,3)
hold on
plot(out.time, out.ta)
plot([0,10],[mins,mins])
plot([0,10],[maxs,maxs])
xlabel('Time (s)')
ylabel('True Air Speed (knots)')
hold off

results.ElevatorAngle = data.elev1__deg(n);
results.IndicatedVelocity = data.Vind_kias(n);
results.TrueVelocity = data.Vtrue_ktas(n);
results.MachRatio = data.Machratio(n);
results.Altitude = data.alt_1ftmsl(n);
results.SinkRate = data.VVI__fpm(n);
results.BankAngle = data.roll__deg(n);
results.FlightPathAngle = data.vpath__deg(n);
end

function filter_cessna(data,name)
    filter_pitch = data.pitch__deg < 5 & data.pitch__deg > -5;
    filter_vy = data.VVI__fpm < 500 & data.VVI__fpm > -500;
    filter_roll = data.roll__deg < 5 & data.roll__deg > -5;
    filter_alt = data.alt_1ftmsl < 5200 & data.alt_1ftmsl > 4800;
    filters = filter_pitch & filter_vy & filter_roll & filter_alt ;

    figure
    hold on
    title(name)
    plot(data.real_time(filters), data.Vtrue_ktas(filters),'.')
    plot([0,data.real_time(end)], [73,73])
    plot([0,data.real_time(end)], [77,77])
    plot([0,data.real_time(end)], [98,98])
    plot([0,data.real_time(end)], [102,102])
    plot([0,data.real_time(end)], [128,128])
    plot([0,data.real_time(end)], [132,132])
    xlabel('Time (s)')
    ylabel('Speed (knots)')
    hold off

end

function filter_747(data,name)
filter_pitch = data.pitch__deg < 5 & data.pitch__deg > -5;
filter_vy = data.VVI__fpm < 500 & data.VVI__fpm > -500;
filter_roll = data.roll__deg < 5 & data.roll__deg > -5;
filter_alt = data.alt_1ftmsl < 20200 & data.alt_1ftmsl > 19800;
filters = filter_pitch & filter_vy & filter_roll & filter_alt ;

figure
hold on
title(name)
plot(data.real_time(filters), data.Vtrue_ktas(filters),'.')
plot([0,data.real_time(end)], [348,348])
plot([0,data.real_time(end)], [352,352])
plot([0,data.real_time(end)], [397,397])
plot([0,data.real_time(end)], [401,401])
plot([0,data.real_time(end)], [448,448])
plot([0,data.real_time(end)], [452,452])
xlabel('Time (s)')
ylabel('Speed (knots)')
hold off
end