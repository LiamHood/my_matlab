close all; clear; clc;

%% Tests
load("Data_ft3.mat")
filters = (Cessna0.vZ__ms)<500 & (Cessna0.vZ__ms)>-500 ...
    & Cessna0.alt_1ftmsl < 5200 & Cessna0.alt_1ftmsl > 4800;
velocity = Cessna0.Vtrue_ktas(filters);
% figure
% plot(Cessna0.Vind_kias,'.')

xvec = [497,578,628,719,787,869,951,1034,1119,1180,1243];

for ii = 1:11
    table2(ii,:) = printItem2(Cessna0,xvec(ii),ii);
end
table2(12,:)  = printItem2(Cessna02,425,12);
disp(table2)




BHP = [51,52,52,53,55,57,62,72,80,90,103,103];
for ii = 1:11
    table3(ii,:) = printItem3(Cessna0,xvec(ii),table2(ii,:),BHP(ii));
end
table3(12,:) = printItem3(Cessna02,425,table2(12,:),BHP(12));

figure
hold on
cdtas = fit(table2(:,5),table2(:,13), 'exp2');
plot(cdtas,table2(:,5),table2(:,13),'*')
cdtas3 = fit(table3(:,3),table3(:,6), 'exp2');
plot(cdtas3,table3(:,3),table3(:,6),'*')
ylabel('CD')
xlabel('TAS')
legend("X-Plane Data","X-Plane Curve Fit", "Calculated Data", "Calculated Curve Fit")
hold off

figure
hold on
claoa = fit(table2(:,7),table2(:,8), 'poly1');
plot(claoa,table2(:,7),table2(:,8),'*')
claoa = fit(table3(:,8),table3(:,7), 'poly1');
plot(claoa,table3(:,8),table3(:,7),'*')
hold off
ylabel('CL')
xlabel('AoA')
legend("X-Plane Data","X-Plane Curve Fit", "Calculated Data", "Calculated Curve Fit")

figure
hold on
cdcl = fit(table2(:,8),table2(:,13), 'exp2');
plot(cdcl,table2(:,8),table2(:,13),'*')
cdcl = fit(table3(:,7),table3(:,6), 'exp2');
plot(cdcl,table3(:,7),table3(:,6),'*')
hold off
ylabel('CD')
xlabel('CL')
legend("X-Plane Data","X-Plane Curve Fit", "Calculated Data", "Calculated Curve Fit")

figure
ldratio = fit(table2(:,5),table2(:,14), 'poly2');
plot(ldratio,table2(:,5),table2(:,14),'*')
xlabel('TAS')
ylabel('L/D')
disp(ldratio)


disp(table3)
function row = printItem2(data,x,num)
    row(1) = num;
    row(2) = data.real_time(x-5);
    row(3) = data.real_time(x);
    row(4) = data.Vind_kias(x);
    row(5) = data.Vtrue_ktas(x);
    row(6) = data.alt_1ftmsl(x);
    row(7) = data.alpha__deg(x);
    row(8) = data.cltotal(x);
    row(9) = data.curnt___lb(x);
    row(10) = data.MP__1_inhg(x);
    row(11) = data.rpm_1_prop(x);
    row(12) = data.thrst_1lb(x);
    row(13) = data.cdtotal(x);
    row(14) = data.LDratio(x);
    fprintf("\nTest Point Number: \t%d \n", num)
    fprintf("Start Time: \t\t%f \n", data.real_time(x-5))
    fprintf("End Time: \t\t%f\n", data.real_time(x))
    fprintf("IAS: \t\t\t%f \n", data.Vind_kias(x))
    fprintf("TAS: \t\t\t%f \n", data.Vtrue_ktas(x))
    fprintf("Altitude: \t\t%f \n", data.alt_1ftmsl(x))
    fprintf("AoA: \t\t\t%f \n", data.alpha__deg(x))
    fprintf("CL: \t\t\t%f \n", data.cltotal(x))
    fprintf("Weight: \t\t%f \n", data.curnt___lb(x))
    fprintf("Manifold Pressure: \t%f \n", data.MP__1_inhg(x))
    fprintf("RPM: \t\t\t%f \n", data.rpm_1_prop(x))
    fprintf("Thrust: \t\t%f \n", data.thrst_1lb(x))
    fprintf("L/D: \t\t\t%f \n", data.LDratio(x))
end

function row = printItem3(data,x,point,BHP)
    qbar = .5*(20.48e-4)*(point(5)*1.68781)^2;
    TAS = point(4)*sqrt(23.77e-4/20.48e-4);
    thrust = 325*.8*BHP/TAS;
    cd = thrust/(qbar*174);
    row(1) = point(1);
    row(2) = data.flaphandl(x);
    row(3) = TAS;
    row(4) = thrust;
    row(5) = qbar;
    row(6) = cd;
    row(7) = point(8);
    row(8) = point(7);
    row(9) = point(8)/cd;
    fprintf("\nTest Point Number: \t%d \n", point(1))
    fprintf("Flap Setting: \t\t%f \n", data.flaphandl(x))
    fprintf("TAS: \t\t\t%f \n", point(4)*sqrt(23.77e-4/20.48e-4))
    fprintf("Drag: \t\t\t%f \n", thrust)
    fprintf("Q_bar: \t\t\t%f \n", qbar)
    fprintf("C_D: \t\t\t%f \n", cd)
    fprintf("C_L: \t\t\t%f \n", data.cltotal(x))
    fprintf("AoA: \t\t\t%f \n", data.alpha__deg(x))
    
end

