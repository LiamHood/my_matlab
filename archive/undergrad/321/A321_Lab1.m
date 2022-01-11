numDataPoints = 30 ;
[t, data] = serial_reader(numDataPoints) ;
disp('Time')
disp(t)
disp( 'Data' )
disp(data)

T0 = 298.15 ;
T_room = 23.6+273.15 ;
B = 3950 ;
R0 = 10000 ;
%T = 1/( (1/T0) + (1/B)*ln(R/R0) )
%% Calibration
resistance = data(:,1) ;
R = mean(resistance) ; %resistance average
som = -(B*T0^2)/(R*(B+T0)*log(R/R0)^2) ; %sensitivity of measurement
temp = data(:,2) ; %temp in celcius
T_meas = mean(temp) ; %average measured temp in C
T_measured = T_meas-273.15 ;
offset = T_meas - T_room ; %Temperature offset
inaccuracy = offset/T_room ; 
accuracy = 1 - inaccuracy ;
T_max = 125 ;
T_min = -55 ;
intervals = 1024 ;
resolution = (T_max-T_min)/intervals ;

%% Linearity
numDataPoints = 30 ;
[t, data] = serial_reader(numDataPoints) ;

figure
plot( t , temp )
xlabel('Time (s)')
ylabel('Temperature (Kelvin)')


range = max(temp)-min(temp) ;
temp_timeconstant = min(temp)+range*.632 ;
ii = 1 ;
    if temp < temp_timeconstant
        ii = ii+1 ;
    else
        timeconstant = t(ii) ;
        disp( 'The time constant is ' )
        disp(timeconstant)
    end



