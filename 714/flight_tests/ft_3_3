close all; clear; clc;

%% Tests
load("Data_ft3.mat")
filters = (Cessna02.vZ__ms)<500 & (Cessna02.vZ__ms)>-500 ...
    & Cessna02.alt_1ftmsl < 5200 & Cessna02.alt_1ftmsl > 4800;
velocity = Cessna02.Vtrue_ktas(filters);
figure
plot(Cessna02.Vind_kias,'.')

% xvec = [461,497,578,628,719,787,869,951,1034,1119,1180,1243];
% 
% for ii = 1:12
%     table2(ii,:) = printItem2(Cessna0,xvec(ii),ii);
% end
% disp(table2)
% 
% 
% 
% 
% BHP = [50,51,52,52,53,55,57,62,72,80,90,103];
% for ii = 1:12
%     table3(ii,:) = printItem3(Cessna0,xvec(ii),table2(ii,:),BHP(ii));
% end
% 
% figure
% hold on
% cdtas = fit(table2(:,13),table2(:,5), 'exp2');
% plot(cdtas,table2(:,13),table2(:,5),'*')
% cdtas3 = fit(table3(:,6),table3(:,3), 'exp2');
% plot(cdtas3,table3(:,6),table3(:,3),'*')
% xlabel('CD')
% ylabel('TAS')
% legend("X-Plane Data","X-Plane Curve Fit", "Calculated Data", "Calculated Curve Fit")
% hold off
% 
% figure
% hold on
% claoa = fit(table2(:,8),table2(:,7), 'poly1');
% plot(claoa,table2(:,8),table2(:,7),'*')
% claoa = fit(table3(:,7),table3(:,8), 'poly1');
% plot(claoa,table3(:,7),table3(:,8),'*')
% hold off
% xlabel('CL')
% ylabel('AoA')
% legend("X-Plane Data","X-Plane Curve Fit", "Calculated Data", "Calculated Curve Fit")
% 
% figure
% hold on
% cdcl = fit(table2(:,13),table2(:,8), 'exp2');
% plot(cdcl,table2(:,13),table2(:,8),'*')
% cdcl = fit(table3(:,6),table3(:,7), 'exp2');
% plot(cdcl,table3(:,6),table3(:,7),'*')
% hold off
% xlabel('CD')
% ylabel('CL')
% legend("X-Plane Data","X-Plane Curve Fit", "Calculated Data", "Calculated Curve Fit")
% 
% 
% disp(table3)
% function row = printItem2(data,x,num)
%     row(1) = num;
%     row(2) = data.real_time(x-5);
%     row(3) = data.real_time(x);
%     row(4) = data.Vind_kias(x);
%     row(5) = data.Vtrue_ktas(x);
%     row(6) = data.alt_1ftmsl(x);
%     row(7) = data.alpha__deg(x);
%     row(8) = data.cltotal(x);
%     row(9) = data.curnt___lb(x);
%     row(10) = data.MP__1_inhg(x);
%     row(11) = data.rpm_1_prop(x);
%     row(12) = data.thrst_1lb(x);
%     row(13) = data.cdtotal(x);
%     fprintf("\nTest Point Number: \t%d \n", num)
%     fprintf("Start Time: \t\t%f \n", data.real_time(x-5))
%     fprintf("End Time: \t\t%f\n", data.real_time(x))
%     fprintf("IAS: \t\t\t%f \n", data.Vind_kias(x))
%     fprintf("TAS: \t\t\t%f \n", data.Vtrue_ktas(x))
%     fprintf("Altitude: \t\t%f \n", data.alt_1ftmsl(x))
%     fprintf("AoA: \t\t\t%f \n", data.alpha__deg(x))
%     fprintf("CL: \t\t\t%f \n", data.cltotal(x))
%     fprintf("Weight: \t\t%f \n", data.curnt___lb(x))
%     fprintf("Manifold Pressure: \t%f \n", data.MP__1_inhg(x))
%     fprintf("RPM: \t\t\t%f \n", data.rpm_1_prop(x))
%     fprintf("Thrust: \t\t%f \n", data.thrst_1lb(x))
% end
% 
% function row = printItem3(data,x,point,BHP)
%     qbar = .5*(20.48e-4)*(point(5)*1.68781)^2;
%     TAS = point(4)*sqrt(23.77e-4/20.48e-4);
%     thrust = 325*.8*BHP/TAS;
%     cd = thrust/(qbar*174);
%     row(1) = point(1);
%     row(2) = data.flaphandl(x);
%     row(3) = TAS;
%     row(4) = thrust;
%     row(5) = qbar;
%     row(6) = cd;
%     row(7) = point(8);
%     row(8) = point(7);
% 
%     fprintf("\nTest Point Number: \t%d \n", point(1))
%     fprintf("Flap Setting: \t\t%f \n", data.flaphandl(x))
%     fprintf("TAS: \t\t\t%f \n", point(4)*sqrt(23.77e-4/20.48e-4))
%     fprintf("Drag: \t\t\t%f \n", thrust)
%     fprintf("Q_bar: \t\t\t%f \n", qbar)
%     fprintf("C_D: \t\t\t%f \n", cd)
%     fprintf("C_L: \t\t\t%f \n", data.cltotal(x))
%     fprintf("AoA: \t\t\t%f \n", data.alpha__deg(x))
% end


