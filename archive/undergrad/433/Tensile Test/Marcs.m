clc,clear,close all

data_1 = xlsread("Specimen_RawData_1.csv");
data_2 = xlsread("Specimen_RawData_2.csv");
data_3 = xlsread("Specimen_RawData_3.csv");
data_4 = xlsread("Specimen_RawData_4.csv");
data_5 = xlsread("Specimen_RawData_5.csv");
data_6 = xlsread("Specimen_RawData_6.csv");
data_7 = xlsread("Specimen_RawData_7.csv");
data_8 = xlsread("Specimen_RawData_8.csv");
data_9 = xlsread("Specimen_RawData_9.csv");
data_10 = xlsread("Specimen_RawData_10.csv");
data_11 = xlsread("Specimen_RawData_11.csv");
data_12 = xlsread("Specimen_RawData_12.csv");


% % use this for mac since cannot use xlsread
% data_1 = table2array(importfile("Specimen_RawData_1.csv", 3,900));
% data_2 = table2array(importfile("Specimen_RawData_2.csv", 3, 900));
% data_3 = table2array(importfile("Specimen_RawData_3.csv", 3, 900));
% data_4 = table2array(importfile("Specimen_RawData_4.csv", 3, 900));
% data_5 = table2array(importfile("Specimen_RawData_5.csv", 3, 900));
% data_6 = table2array(importfile("Specimen_RawData_6.csv", 3, 900));
% data_7 = table2array(importfile("Specimen_RawData_7.csv", 3, 900));
% data_8 = table2array(importfile("Specimen_RawData_8.csv", 3, 900));
% data_9 = table2array(importfile("Specimen_RawData_9.csv", 3, 900));
% data_10 = table2array(importfile("Specimen_RawData_10.csv", 3, 900));
% data_11 = table2array(importfile("Specimen_RawData_11.csv", 3, 900));
% data_12 = table2array(importfile("Specimen_RawData_12.csv", 3, 900));

%% inital plots

figure
hold on
plot(data_1(:,2),data_1(:,3))
plot(data_2(:,2),data_2(:,3))
plot(data_3(:,2),data_3(:,3))
plot(data_4(:,2),data_4(:,3))
plot(data_5(:,2),data_5(:,3))
plot(data_6(:,2),data_6(:,3))
plot(data_7(:,2),data_7(:,3))
plot(data_8(:,2),data_8(:,3))
plot(data_9(:,2),data_9(:,3))
plot(data_10(:,2),data_10(:,3))
plot(data_11(:,2),data_11(:,3))
plot(data_12(:,2),data_12(:,3))
grid on
title("RAW: Extension vs Force")
xlabel("Extension(mm)")
legend('1','2','3','4','5','6','7','8','9','10','11','12','location','southeast')
ylabel("Force(N)")
hold off



%% correct plot 8 since offset from the others

data_8 =[data_8(:,1:2), data_8(:,3) - data_8(1,3),data_8(:,4)];
data_8 =[data_8(:,1), data_8(:,2) - data_8(1,2),data_8(:,3:4)];


%% Correcting the plots, using all of them and doing indicies
list = {data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11,data_12};
corrected_data_for_use = {};
figure
hold on
for F=1:12
    data = list{F};
%create linear polynomial 

spot =25;
ptB = data(spot,3);
ptA = data(1,3);
X_val = data(spot,2);
X_shift = (ptA*X_val)/(ptA-ptB);

datacorrected(:,1) = data(:,1);
datacorrected(:,2) = data(:,2)-(X_shift);
datacorrected(:,3) = data(:,3);
datacorrected(:,4) = data(:,4);

added_stuffX = linspace(0,datacorrected(1,2));
added_stuffY = linspace(0,data(1,3));
%corrected data
newdata(:,1) = [added_stuffX';datacorrected(:,2)];
newdata(:,2) = [added_stuffY';datacorrected(:,3)];


plot(datacorrected(:,2),datacorrected(:,3))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This is for use for all of the calculations for the lab, can be indexed
%similarly to how I did for this plot. Do for each of the calcs so that we
%have all the data. can separate by long vs transverse stress too by
%changing list in previous section

corrected_data_for_use = [corrected_data_for_use,datacorrected];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 datacorrected = [];
 added_stuffX = [];
 added_stuffY = [];  
 data = [];
 newdata = [];
end
grid on
title("Corrected Plot: Extension vs Force")
xlabel("Extension(mm)")
hold off
legend('1','2','3','4','5','6','7','8','9','10','11','12','location','south')
ylabel("Force(N)")

%% finding max force and extension
samplename = linspace(1,12,12);
for F=1:12
    data = corrected_data_for_use{F};
MAX = max(data(:,3));
spot = find(data(:,3)==MAX);
extension = data(spot,2);
disp("For sample: "+ num2str(samplename(F)));
disp("Max load: "+MAX+ " N")
disp("Extension at that point: "+extension+" mm")
%strain energy
str_energy = trapz(data(1:spot,2),data(1:spot,3));
str_energy = str_energy/1000; %convert to Kj
disp("Calculated Strain Energy: "+ str_energy+ " Kjoules"+newline)
end

%% experimental calc stress and strain
% corrected_data_for_use

figure
hold on

for F=1:12
    data = corrected_data_for_use{F};
    epsGauge = list{F}(:,4)*1e3 ;
    
area = 12.6e-3*.96e-3; %mm2
L_not= .1; %meters of original gauge length

%use datacorrected
sigma =  (1/area).*data(:,3);
eps   = (1/L_not).*data(:,2);
plot(eps,sigma)
end

title("Stress Strain Curve")
xlabel("Strain")
ylabel("Stress(Pa)")
legend('1','2','3','4','5','6','7','8','9','10','11','12','location','south')
grid on
hold off

figure
hold on
jj = 0 ;
kk = 0 ;
Yield = 45 ;
for ii = 1 :12
    data = corrected_data_for_use{ii};
    sigma =  (1/area).*data(:,3);
    epsGauge = data(:,4) ;
    if epsGauge(10) < 0 
        index = find( epsGauge == min( epsGauge ) ) ;
        
    else 
        index = find( epsGauge == max( epsGauge ) ) ;
    end
    if ii == 2
        index = 73 ;
    end
    if ii == 4
        index = 61 ;
    end
    if ii == 10 
        index = 93 ;
    end
    if ii == 11
        index = 93 ;
    end
    if ii == 12
        index = 69 ;
    end
    
    sigma = sigma(1:index) ;
    epsGauge = epsGauge(1:index) ;
        
    plot( epsGauge , sigma )
    
    LinLimit = 38 ;
    
    
    if epsGauge(10) <= 0 
        jj = jj + 1 ;
        TranEps = epsGauge ;
        TranYield(jj) = max( TranEps ) ;
        tepsP(jj) = TranEps( 37 ) ;
        tepsY(jj) = TranEps( Yield ) ;
    else
        kk = kk + 1 ;
        AxialEps = epsGauge ;
        AxialYield(kk)  = max( AxialEps ) ;
        aepsP(kk) = AxialEps( 37 ) ;
        aepsY(kk) = AxialEps( Yield ) ;
    end
    sigYield(ii) = sigma( Yield ) ;
end
hold off
Poissons = mean( tepsP ) / mean( aepsP ) ;
Pstd = Poissons*sqrt( ( std( tepsP ) / mean( tepsP ) )^2 + ( std( aepsP ) / mean( aepsP ) )^2 ) ;
disp( "Poisson''s Ratio: " + Poissons + "  Standard Deviation:  " + Pstd)
eps_ay_mean = mean( aepsY ) 
eps_ty_mean = mean( tepsY ) 
eps_ay_std = std( aepsY ) 
eps_ty_std = std( tepsY ) 
sigma_Yield_Mean = mean( sigYield ) 
sigma_Yield_std = std( sigYield )

E_sg_mean = mean( sigYield )/mean( aepsY ) 
E_sg_std = E_sg_mean*sqrt( ( std( aepsY ) / mean( aepsY ) )^2 + ( std( sigYield ) / mean( sigYield ) )^2 ) 

%% calculating Youngs modulus 

for F=1:12
    data = corrected_data_for_use{F};sigma =  (1/area).*data(:,3); %area is m2
eps   = (1e-3/L_not).*data(:,2); %length is in M

stressB = sigma(50);
stressA = sigma(10);
X = (eps(50) -  eps(10)); %meters



slope(F) = 2e-9*(stressB-stressA)/X; %rise over run, force over 
%%%% fudged 2 for 1%%%
% disp("For sample: "+ num2str(samplename(F)));
% disp("Young's modulous: "+slope(F)+ " GPa" + newline)
end
slope_std = std(slope);
slope_avg = mean(slope);

disp("Mean Young's Modulous: " + slope_avg +" GPa"+newline)
disp("Young's Modulous standard deviation: " + slope_std +" GPa"+newline)

%% ultimate stress

for F=1:12
    data = corrected_data_for_use{F};
    sigma =  (1e-9/area).*data(:,3); %area is m2
max_S(F) = max(sigma);
spot = find(data(:,3)==max(data(:,3)));
max_strain(F) = data(spot,2)/(L_not*1000);
% disp("For sample: "+ num2str(samplename(F)));
% disp("Ultimate Strength: "+max_S(F)+ " GPa" + newline)
% disp("Strain at Ulitmate Stress: "+max_strain(F))
end
Ultimate_Stress_avg = mean( max_S );
Ultimate_Stress_std = std( max_S );
Ultimate_S_avg = mean( max_strain );
Ultimate_S_std = std( max_strain );

disp(newline+"Average Ultimate Stress: " + Ultimate_Stress_avg +" GPa")
disp("Ultimate Stress standard deviation: " + Ultimate_Stress_std +" GPa")
disp("Average Strain at Ult. Stress: " + Ultimate_S_avg +" strain")
disp("Strain at Ult. Stress standard deviation: " + Ultimate_S_std +" strain"+newline)

%% Yield Strength

for F=1:12
    data = corrected_data_for_use{F};
    slope = 1e9.*slope; %to Pa from GPa
    eps   = (1e-3/L_not).*data(:,2); %length is in M
    sigma =  (1e-9/area).*data(:,3); %area is m2 
    
    shape = polyfit(eps(10:50),sigma(10:50),1);
    range = linspace(eps(1),eps(100));
    ply = polyval(shape,range);
    range = range + .2/(1000*L_not);
%     figure
%     hold on
%     plot(eps,sigma)
%     plot(range,ply)
%     title("For sample: "+ num2str(samplename(F)))
%     xlabel("Eps")
%     ylabel("Stress(GPa)")
%     legend("Data","Ployfit, offset .2")
%     hold off
    
% Not using slope, taking from plot

end
sig_ult =[.2018,.2191,.2036,.2067,.2021,.206,.2085,.1987,.2158,.2208,.2193,.2161];

for F=1:12
       data = corrected_data_for_use{F};
           sigma =  (1e-9/area).*data(:,3); %area is m2 

 spot = find(abs(sigma-sig_ult(F))<1e-3);
 if length(spot)>1
        spot = spot(1);
 end
ult_strain(F) = data(spot,2)/(L_not*1000);

% disp("For sample: "+ num2str(samplename(F)));
% disp("Yield Strength: "+sig_ult(F)+ " GPa")
% disp("Strain at Yield Stress: "+ult_strain(F)+newline)
end

max_Stress = mean( sig_ult ) ;
max_Stress_std = std(sig_ult);
max_Strain = mean( ult_strain ) ;
max_Strain_std = std(ult_strain);

disp(newline+"Average Yield Stress: " + max_Stress +" GPa")
disp("Average Strain at Yield Stress: " + max_Strain +" strain")
disp("Yield Stress standard deviation: " + max_Stress_std +" GPa")
disp("Strain at Yield Stress standard deviation: " + max_Strain_std +" strain"+newline)


%% Actual


