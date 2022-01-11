clear ; close all ; clc ;
length = 215 ;
width = 61 ;
thick = 1.7 ;
for ii = 1:15
    rawname = [ 'Specimen_RawData_' , num2str( ii ) , '.csv' ] ;
    data{ii} = xlsread( rawname ) ;
    loadFull = data{ii}(:,3) ;
    load(ii) = 10 ;
    interest = find( loadFull > 9.95 & loadFull < 10.5 ) ;
    extension(ii) = mean( data{ii}(interest,2) ) ;
    sg(ii) = mean( data{ii}(interest,4) ) ;
end
stressm = sg.*69e3 ;
for ii = 1:4
    dfc(ii) = 215 - 72 + 3 + 3*15 - (ii-1)*15 ;
end
for ii = 5:9
    dfc(ii) = 215 - 72 - 19 - 3 - (ii-5)*15 ;
end
dfc(10) = 2.8 ;
dfca = [ dfc , 210 ] ;
I = (1/12)*(width*1e-3)*(thick*1e-3)^3 ;
c = (thick*1e-3)/2 ;
M = 10.*(.215 - dfca*1e-3 ) ;
stresst = -(M.*c./I)*1e-6 ;

n = 2 ;
p = polyfit( dfc , sg(1:10) , n) ;
sgfit = 0 ;
dfcfit = linspace( 0 , 210 , 1e2 ) ;
for ii = 1:n+1
    sgfit = sgfit + p(ii)*dfcfit.^(n-(ii-1)) ;
end

figure
hold on
plot( dfc , sg(1:10) , '*' )
plot( dfcfit , sgfit )
hold off
xlabel( 'Distance from Clamp [mm]' )
ylabel( 'Strain' )
title( 'Measured Strain along Damaged Beam' )
legend( 'Strain Gauge Values' , 'Best Fit Curve' , 'Location' , 'southeast' )

figure
hold on
plot( dfc , stressm(1:10) , '*' )
plot( dfcfit , sgfit*69e3 )
plot( dfca , stresst )
hold off
xlabel( 'Distance from Clamp [mm]' )
ylabel( 'Stress [MPa]' )
title( 'Stress in the Beam' )
legend( 'Stress from Strain Gauge Values' , 'Best Fit Curve' , 'Stress for an Undamaged Beam' , 'Location' , 'southeast' )

for ii = 1:5
    dfh(ii) = ( 61 - 19 )/2 - ( 2.4 + 3.9*(ii-1) ) ;
end
sgtofit = [ sg(11:13) , sg( 15 ) ] ;
dfhtofit = [ dfh(1:3) , dfh( 5 ) ] ;
p = polyfit( dfhtofit , sgtofit , n) ;
sgsfit = 0 ;
dfhfit = linspace( 0 , 20 , 1e2 ) ;
for ii = 1:n+1
    sgsfit = sgsfit + p(ii)*dfhfit.^(n-(ii-1)) ;
end
figure
hold on
plot( dfh , sg(11:15) , '*' )
plot( dfhfit , sgsfit )
hold off
xlabel( 'Distance from Damage [mm]' )
ylabel( 'Strain' )
title( 'Measured Strain Next to Damaged' )
legend( 'Strain Gauge Values' , 'Best Fit Curve' , 'Location' , 'southeast' )

figure
hold on
plot( dfh , sg(11:15)*69e3 , '*' )
plot( dfhfit , sgsfit*69e3 )
hold off
xlabel( 'Distance from Damage [mm]' )
ylabel( 'Stress [MPa]' )
title( 'Stress from Strain Gauges Next to Damaged' )
legend( 'Strain Gauge Values' , 'Best Fit Curve' , 'Location' , 'southeast' )
