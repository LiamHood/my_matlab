%% Lab 2
% Aero 321
% Liam Hood
clear; clc; close all;
%% Slow
clear
load('distance.mat')
load('time.mat')
distance = data(:,1) ;

distance_avg = mean(distance) ;

dev = distance - mean(distance) ;
var_top = 0 ;
for ii = 1:length(distance)
    var_top = dev(ii)^2 + var_top ;
end
var = var_top / ( length(distance) - 1 ) ;
sigma = sqrt(var) ; %standard deviation

alpha = sigma/(sqrt(length(distance))); %Standard error
distance_actual = 10 ; %cm
accuracy = 1-(distance_actual-distance_avg)/distance_actual ; 
precision = mean( abs(dev) ); 
[ h,p,stats ] = chi2gof( distance ) ;

figure
plot( t , distance )
xlabel('Time (s)')
ylabel('Distance (cm)')
title('Raw Data')

%% Fast
clear
load('data_100Hz.mat')
distance_100 = data_100(:,1) ;

distance_avg_100 = mean(distance_100) ;

dev_100 = distance_100 - mean(distance_100) ;
var_top_100 = 0 ;
for ii = 1:length(distance_100)
    var_top_100 = dev_100(ii)^2 + var_top_100 ;
end
var_100 = var_top_100 / ( length(distance_100) - 1 ) ;
sigma_100 = sqrt(var_100) ; %standard deviation

alpha_100 = sigma_100/(sqrt(length(distance_100))); %Standard error
distance_actual_100 = 10.5 ; %cm
accuracy_100 = 1-(distance_actual_100-distance_avg_100)/distance_actual_100 ; 
precision_100 = mean( abs(dev_100) ); 
[ h,p,stats_100 ] = chi2gof( distance_100 ) ;

figure
plot( t_100 , distance_100 )
xlabel('Time (s)')
ylabel('Distance (cm)')
title('Raw Data at 100 Hz')

%% Filtered
clear
% alpha = .5
load( 'Filtered_alpha_.5.mat' )
distance = data(:,1) ;
distance_f = data(:,2) ;

distance_avg = mean(distance) ;
distance_avg_f = mean(distance_f) ;

dev = distance - mean(distance) ;
dev_f = distance_f - mean(distance_f) ;
var_top = 0 ;
var_top_f = 0 ;
for ii = 1:length(distance)
    var_top = dev(ii)^2 + var_top ;
end
for ii = 1:length(distance_f)
    var_top_f = dev_f(ii)^2 + var_top_f ;
end
var = var_top / ( length(distance) - 1 ) ;
sigma = sqrt(var) ; %standard deviation
var_f = var_top_f / ( length(distance_f) - 1 ) ;
sigma_f = sqrt(var_f) ; %standard deviation

alpha = sigma/(sqrt(length(distance))); %Standard error
distance_actual = 23.5 ; %cm
accuracy = 1-(distance_actual-distance_avg)/distance_actual ; 
precision = mean( abs(dev) ); 
[ h,p,stats ] = chi2gof( distance ) ;

alpha_f = sigma_f/(sqrt(length(distance_f))); %Standard error
distance_actual = 23.5 ; %cm
accuracy_f = 1-(distance_actual-distance_avg_f)/distance_actual ; 
precision_f = mean( abs(dev_f) ); 
[ h,p,stats_f ] = chi2gof( distance_f ) ;

figure
plot( t , distance , t , distance_f)
xlabel('Time (s)')
ylabel('Distance (cm)')
title('Data')
legend( 'Raw' , 'Filtered Alpha = .5' )

%% Moving

[t_move, data_move] = serial_reader(100) ;

figure
plot( t_move , data_move(:,1) , t , data_move(:,2))
xlabel('Time (s)')
ylabel('Distance (cm)')
title('Changing Distance Data')
legend( 'Raw' , 'Filtered Alpha = .5' )

%% Moving many alphas

[t_filter, data_filter] = serial_reader(100) ;

figure
for ii = 1:6
    hold on 
    plot( t_filter , data_filter(:,ii) )
end
xlabel('Time (s)')
ylabel('Distance (cm)')
title('Changing Distance Many Filters Data')
legend( 'Raw' , 'Filtered Alpha = .1' , 'Filtered Alpha = .3' , 'Filtered Alpha = .5' , 'Filtered Alpha = .7' , 'Filtered Alpha = .9')

%% Proof Square
clear
%[t_l, data_l] = serial_reader(100) ;
%[t_w, data_w] = serial_reader(100) ;
load( 'ProofSquareL.mat' )
load( 'ProofSquareW.mat' )
% error in Length calc
    deviance_l = zeros( length( data_l ) , 6 ) ;
    variance_l = zeros( 1 ,  6 ) ;
    sigma_l = zeros( 1 , 6 ) ;
    for jj = 1:6
        deviance_l(:,jj) = data_l(:,jj) - ones(length(data_l),1)*mean( data_l(:,jj) ) ;
        variance_l(jj) = sum( deviance_l(:,jj).^2 )/(length(deviance_l(:,jj)) - 1 ) ;
        sigma_l(jj) = sqrt( variance_l(jj) ) ;
    end
err_l = sigma_l*1.96 ;

% error in Width calc
    deviance_w = zeros( length( data_w ) , 6 ) ;
    variance_w = zeros(1 , 6 ) ;
    sigma_w = zeros( 1 , 6 ) ;
    for jj = 1:6
        deviance_w(:,jj) = data_w(:,jj) - ones(length(data_w),1)*mean( data_w(:,jj) ) ;
        variance_w(jj) = sum( deviance_w(:,jj).^2 )/(length(deviance_w(:,jj)) - 1 ) ;
        sigma_w(jj) = sqrt( variance_w(jj) ) ;
    end
err_w = sigma_w*1.96 ;

for kk = 1:6
    error(kk) = sqrt( err_w(kk)^2 + err_l(kk)^2 );
    area(kk) = mean( data_w(:,kk) ) * mean( data_l(:,kk) ) ;
end

x = [ 1 .1 .3 .5 .7 .9 ] ;

figure
hold on
    plot( x , area , '.' )
    errorbar( x , area , error , 'o')
hold off
xlabel('Time (s)')
ylabel('Area (cm)')
title('Proof Square Measurements')

load( 'ProofSquareL.mat' )
load( 'ProofSquareW.mat' )

% Square 2
load( 'ProofSquare2L.mat' )
load( 'ProofSquare2W.mat' )

% error in Length calc
    deviance2_l = zeros( length( data2_l ) , 6 ) ;
    variance2_l = zeros( 1 ,  6 ) ;
    sigma2_l = zeros( 1 , 6 ) ;
    for jj = 1:6
        deviance2_l(:,jj) = data2_l(:,jj) - ones(length(data2_l),1)*mean( data2_l(:,jj) ) ;
        variance2_l(jj) = sum( deviance2_l(:,jj).^2 )/(length(deviance2_l(:,jj)) - 1 ) ;
        sigma2_l(jj) = sqrt( variance2_l(jj) ) ;
    end
err2_l = sigma2_l*1.96 ;

% error in Width calc
    deviance2_w = zeros( length( data2_w ) , 6 ) ;
    variance2_w = zeros(1 , 6 ) ;
    sigma2_w = zeros( 1 , 6 ) ;
    for jj = 1:6
        deviance2_w(:,jj) = data2_w(:,jj) - ones(length(data2_w),1)*mean( data2_w(:,jj) ) ;
        variance2_w(jj) = sum( deviance2_w(:,jj).^2 )/(length(deviance2_w(:,jj)) - 1 ) ;
        sigma2_w(jj) = sqrt( variance2_w(jj) ) ;
    end
err2_w = sigma2_w*1.96 ;

for kk = 1:6
    error2(kk) = sqrt( err2_w(kk)^2 + err2_l(kk)^2 );
    area2(kk) = mean( data2_w(:,kk) ) * mean( data2_l(:,kk) ) ;
end

x = [ 1 .1 .3 .5 .7 .9 ] ;

figure
hold on
    plot( x , area2 , '.' )
    errorbar( x , area2 , error2 , 'o')
hold off
xlabel('Time (s)')
ylabel('Area (cm)')
title('Proof Square 2 Measurements')
%% 2
load( 'ProofSquare2L.mat' )
load( 'ProofSquare2W.mat' )
area2 = data2_l .* data2_w ;
    deviance2 = zeros( length( area2 ) , 6 ) ;
    variance2 = zeros( length( area2 ) , 6 ) ;
    sigma2 = zeros( length( area2 ) , 6 ) ;
    for jj = 1:6
        deviance2(:,jj) = area2(:,jj) - ones(length(area2),1)*mean( area2(:,jj) ) ;
        variance2(:,jj) = sum( deviance2(:,jj).^2 )/(length(deviance2(:,jj)) - 1 ) ;
        sigma2(:,jj) = sqrt( variance2(:,jj) ) ;
    end
err2 = sigma2*1.96 ;

figure
for ii = 1:6
    hold on 
    plot( t2_w , area2(:,ii) )
    errorbar( t2_w , area2(:,ii) , err2(:,ii) )
end
xlabel('Time (s)')
ylabel('Area (cm)')
title('Second Proof Square Measurements')
legend( 'Raw' , 'Filtered Alpha = .1' , 'Filtered Alpha = .3' , 'Filtered Alpha = .5' , 'Filtered Alpha = .7' , 'Filtered Alpha = .9')



