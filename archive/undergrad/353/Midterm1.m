%% Midterm 1
% Spacecraft Environments
% Liam Hood
function Midterm1()
clc; clear; close all;

AtmosphericTable = [ 80	0.021999025	0.024571135	3.22009E-20	4.11549E-20 ; ...
100	0.021546553	0.024570958	3.06282E-20	4.1155E-20 ; ...
120	0.020572841	0.024107455	2.72417E-20	3.95418E-20 ; ...
140	0.020151824	0.024107237	2.57837E-20	3.95419E-20 ; ...
160	0.02086041	0.025089884	2.82738E-20	4.29638E-20 ; ...
180	0.020421949	0.025089881	2.67765E-20	4.29638E-20 ; ...
200	0.022634086	0.026635973	3.4598E-20	4.83512E-20 ; ...
220	0.021613262	0.025956358	3.11248E-20	4.59902E-20 ; ...
240	0.021135058	0.025173841	2.95466E-20	4.32741E-20 ; ...
260	0.021875748	0.025149913	3.2131E-20	4.31986E-20 ; ...
280	0.020375056	0.024860903	2.70977E-20	4.22018E-20; ...
300	0.018972346	0.024833729	2.24304E-20	4.21212E-20; ...
320	0.018104496	0.02542809	1.97679E-20	4.42149E-20; ...
340	0.017154812	0.025014337	1.70363E-20	4.27978E-20; ...
360	0.016438595	0.024711836	1.60011E-20	4.17767E-20; ...
380	0.01553743	0.024121928	1.40283E-20	3.97589E-20; ...
400	0.014143682	0.02359261	1.3322E-20	3.79628E-20; ...
420	0.01066863	0.022535034	1.24689E-20	3.43184E-20; ...
440	0.007258539	0.021449041	1.06353E-20	3.05727E-20; ...
460	0.00432235	0.020404336	9.03802E-21	2.69659E-20; ...
480	0.00351503	0.019459819	8.28983E-21	2.37016E-20; ...
500	0.003010642	0.018651911	7.961E-21	2.09065E-20; ...
520	0.002177618	0.018116948	7.74687E-21	1.90794E-20; ...
540	0.001817642	0.017637395	7.79726E-21	1.74364E-20; ...
560	0.001737794	0.017266981	8.04958E-21	1.61774E-20; ...
580	0.001645458	0.017042738	8.47109E-21	1.54638E-20; ...
600	0.001588727	0.016836129	9.20218E-21	1.48111E-20 ] ;

AODensity = [ 200	1.00E+16	1.00E+16 ; ... 
220	6.00E+15	8.00E+15;...
240	3.00E+15	6.00E+15;...
260	1.00E+15	4.00E+15;...
280	7.00E+14	3.00E+15;...
300	5.00E+14	2.00E+15;...
320	3.00E+14	1.00E+15;...
340	2.00E+14	8.00E+14;...
360	1.00E+14	6.00E+14;...
380	7.00E+13	5.00E+14;...
400	3.00E+13	4.00E+14;...
420	8.00E+12	3.80E+14;...
440	3.00E+12	3.60E+14;...
460	1.00E+12	3.40E+14;...
480	7.00E+11	3.20E+14;...
500	5.00E+11	3.00E+14;...
520	2.00E+11	2.60E+14;...
540	9.00E+10	2.30E+14;...
560	6.00E+10	2.00E+14;...
580	3.00E+10	1.60E+14;...
600	1.00E+10	1.30E+14 ] ;

g = 9.81 ; % m*s^-2
R = 287 ; % J*mol^-1*K^-1
kB = 138e-23 ; % J/K
Na = 6.0221409e23 ; % avogadro's number
mu = 398600 ;
re = 6378 ;

h = AtmosphericTable( 1:17 , 1 ) ;
Mmin = AtmosphericTable( 1:17 , 2 ) ;
Mmax = AtmosphericTable( 1:17 , 3 ) ;
Mavg = ( Mmin + Mmax ) ./ 2 ;

%% 1
sminT = linspace( 180 , 700 , 17 ) ;
smaxT = linspace( 180 , 1800 , 17 ) ;
savgT = ( sminT + smaxT ) ./ 2 ;

psmin = 1 ; % [Pa] at 80 km
psavg = 1 ; 
psmax = 1 ; % [Pa] at 80 km
for ii = 2:17
    psmin(ii) = psmin(ii-1)*exp( -(Mmin(ii)*g*20e3) / (R*sminT(ii)) ) ;
    psavg(ii) = psavg(ii-1)*exp( -(Mavg(ii)*g*20e3) / (R*savgT(ii)) ) ;
    psmax(ii) = psmax(ii-1)*exp( -(Mmax(ii)*g*20e3) / (R*smaxT(ii)) ) ;
end
figure
plot( psmin , h , psmax , h , psavg , h) 
title( 'Pressure in the Upper Atmosphere' ) 
xlabel( 'Pressure (Pa)' )
ylabel( 'Altitude (km)' )
legend( 'Solar Minimum' , 'Solar Maximum' , 'Solar Average' )
% The pressure drops less with altitude for solar maximum because the
% molecules are hotter so moving more, creating more pressure expanding the
% atmosphere

% particle density
ng_smin = psmin./(kB.*sminT) ;
ng_savg = psavg./(kB.*savgT) ;
ng_smax = psmax./(kB.*smaxT) ;
figure
plot( ng_smin , h , ng_smax , h , ng_savg , h ) 
title( 'Particle Density in the Upper Atmosphere' ) 
xlabel( 'Particle Density (particle/m^3)' )
ylabel( 'Altitude (km)' )
legend( 'Solar Minimum' , 'Solar Maximum' , 'Solar Average' )
% The particle density decreases with altitude. It is faster for solar
% maximum because the molecules are warmer and collide more often expanding
% the space between molecules

% Thermal velocity
for ii = 1:17
    v_smin(ii) = sqrt( 8*kB*sminT(ii)/(pi*(Mmin(ii)/Na)) ) ;
    v_savg(ii) = sqrt( 8*kB*savgT(ii)/(pi*(Mavg(ii)/Na)) ) ;
    v_smax(ii) = sqrt( 8*kB*smaxT(ii)/(pi*(Mmax(ii)/Na)) ) ;
end
figure
plot( v_smin , h , v_smax , h , v_savg , h ) 
title( 'Thermal Velocity in the Upper Atmosphere' ) 
xlabel( 'Thermal Velocity (m/s)' )
ylabel( 'Altitude (km)' )
legend( 'Solar Minimum' , 'Solar Maximum' , 'Solar Average' )
% The thermal velocity increases with temperature which makes sense. And
% they start at about the same velocity because the temperature and
% composition of the atmosphere is not very effected at lower altitudes. At
% the high altitude the solar maximum thermal velocity is much higher
% because the molecules are much warmer. It is weird that the curve is less
% smooth

% Mean free path
for ii = 1:17
    mfp_smin = 1./(4.*ng_smin*(AtmosphericTable(ii,4))) ;
    mfp_savg = 1./(4.*ng_savg*((AtmosphericTable(ii,4)+AtmosphericTable(ii,5))./2)) ;
    mfp_smax = 1./(4.*ng_smax*(AtmosphericTable(ii,5))) ;
end
figure
plot( mfp_smin , h , mfp_smax , h , mfp_savg , h ) 
title( 'Mean Free Path in the Upper Atmosphere' ) 
xlabel( 'Mean Free Path' )
ylabel( 'Altitude (km)' )
legend( 'Solar Minimum' , 'Solar Maximum' , 'Solar Average' )
% It makes sense that the mean free path increase with altitude since the
% particel density of the atmosphere decreases. It also makes sense that
% solar maximum makes it shorter because the atmosphere expands more with
% the extra heat.

%% 2
% p = rho*R*T
% a 
n = 1e2 ;
p0 = 1 ; % avg pressure [Pa] at 80km
T0 = 273.15 - 80 ; %avg temp [K] at 80 km
L = 3 ; % [K/km]
p = p0.*( 1 + (L/T0).*h ).^((Mavg.*g)./(R.*L)) ;
T = T0 + L/h ;
for ii = 1:length(h)
    rho(ii) = p(ii)/(R*T(ii)) ;
end
figure
plot( rho , h ) 
title( 'Density in the Upper Atmosphere' ) 
xlabel( 'Density (kg/m^3)' )
ylabel( 'Altitude (km)' )

% b
a = 400 + re ;
m = 80 ;
A = 1 ;
Cd = 2.2 ;
time = 0 ;
ii = 1 ;
while a(ii) > ( 80 + re )
    adot(ii) = -((rho(ii)*A*Cd)/m)*((mu*a(ii))^.5)*(1e3) ;
    Life(ii) = -20/adot(ii) ;
    time(ii+1) = sum(Life) ;
    a(ii+1) = a(ii) - 20 ;
    ii = ii + 1 ;
end
disp([ 'The time to de-orbit is ' , num2str(time(ii)) ])

% c
m = 80 ;
A = linspace( 1 , 100 , 100 );
Cd = 2.2 ;
time = 0 ;
for jj = 1:100 
    a = 400 + re ;
    ii = 1 ;
    while a(ii) > ( 80 + re )
        adot(ii,jj) = -((rho(ii)*A(jj)*Cd)/m)*((mu*a(ii))^.5)*(1e3) ;
        Life(ii) = -20/adot(ii,jj) ;
        time(jj) = sum(Life) ;
        a(ii+1) = a(ii) - 20 ;
        ii = ii + 1 ;
    end
end
figure
plot( A , time )
title( 'S/c lifespan with Different Sails' )
xlabel( 'Sail Area (m^2)' )
ylabel( 'Lifespan' )

%% 3 
mu = 398600 ;
re = 6378 ;
v = sqrt( mu / re )*10^5 ; % cm/s
E_ag = 10.5e-24 ; % cm^3 from slides
phi_min = AODensity(12,2)*10^-6*v ;
phi_max = AODensity(12,3)*10^-6*v ;
erate_min = E_ag*phi_min ;
erate_max = E_ag*phi_max ;
edepth_min = 0 ;
t = 0 ;
ii = 2 ;
while edepth_min(ii-1) <= .1 
    edepth_min(ii) = edepth_min(ii - 1) + erate_min*10^5 ; % erosion since start
    t(ii) = t(ii-1)+10^5 ; % seconds since start
    ii = ii + 1 ;
end
edepth_max = 0 ;
tmax = 0 ;
ii = 2 ;
while edepth_max(ii-1) <= .1 
    edepth_max(ii) = edepth_max(ii - 1) + erate_max*10^5 ; % erosion since start
    tmax(ii) = tmax(ii-1)+10^5 ; % seconds since start
    ii = ii + 1 ;
end

figure
plot( t/(60*60*24) , edepth_min , tmax/(60*60*24) , edepth_max )
title( 'Erosion of Silver' ) 
xlabel( 'Time (days)' )
ylabel( 'Erosion Depth (cm)' )
legend( 'Solar Minimum' , 'Solar Maximum' )

%% 4 
TML = .91 ;
CVCM = .07 ;
% a 
temp = linspace( 50 , 300 , 1e3 ) ;
Ea = [ 4 , 8 , 12 , 16 ] ;
tau0 = 5e-12 ;
for ii = 1:length(Ea) 
    tau(ii,:) = tau0.*exp(Ea(ii)./(R.*temp)) ;
end
plot( temp , tau(1,:) , temp , tau(2,:) , temp , tau(3,:) , temp , tau(4,:) )
title( 'Residence Time vs Temperature at Different Activation Energies' ) 
xlabel( 'Temperature (Kelvin)' )
ylabel( 'Residence Time (s)' )
legend( 'Ea = 4' , 'Ea = 8' , 'Ea = 12' , 'Ea = 16' )

% b
m0 = 11 ; % grams
t = 24 ;
temp = 273.15 + [ 0 , 40 , 80 ] ; 
Ea = 12*4184 ; 
q0 = ( m0*TML )./(2*exp(-Ea./(R*398*sqrt(24)))) ;
for jj = 1:length(temp)
    m0 = 11 ; % grams
    for ii = 1:7
        t(ii) = (ii) ;
        q0 = ( m0*TML )./(2*exp(-Ea./(R*398*sqrt(24)))) ;
        deltam(ii,jj) = 2.*q0.*exp(-Ea./(R.*temp(jj))).*(sqrt(t(ii))-sqrt(t(ii)-1)) ;
        m0 = m0 - deltam(ii,jj) ;
    end
end
day = t ;
figure 
plot( day , deltam(:,1) , day , deltam(:,2) , day , deltam(:,3) )
title( 'Daily Mass Loss' )
xlabel( 'Mission Day' )
ylabel( 'Mass Loss (g)' )
legend( '0 degrees C' , '40 degrees C' , '80 degrees C' )

% c
F = linspace( .00033 , .044 , 5 ) ;
temp = 274 ;
Ea = 10*4184 ; 
m0 = 11 ;
for jj = 1:length(F)
    m0 = 11 ; % grams
    for ii = 1:7
        t(ii) = (ii) ;
        q0 = ( m0*TML )./(2*exp(-Ea./(R*398*sqrt(24)))) ;
        deltam(ii,jj) = 2.*q0.*exp(-Ea./(R.*temp)).*(sqrt(t(ii))-sqrt(t(ii)-1)) ;
        deltamd(ii,jj) = F(jj)*deltam(ii,jj) ;
        m0 = m0 - deltam(ii,jj) ;
    end 
end
day = t ;
figure 
plot( day , deltam(:,1) , day , deltam(:,2) , day , deltam(:,3) , day , deltam(:,4) , day , deltam(:,5) )
title( 'Daily Mass Loss' )
xlabel( 'Mission Day' )
ylabel( 'Mass Loss (g)' )
legend( [ 'View Factor = ' , num2str(F(1)) ] , [ 'View Factor = ' , num2str(F(2)) ] , [ 'View Factor = ' , num2str(F(3)) ] , [ 'View Factor = ' , num2str(F(4)) ] , [ 'View Factor = ' , num2str(F(5)) ] ) ;
figure 
plot( day , deltamd(:,1) , day , deltamd(:,2) , day , deltamd(:,3) , day , deltamd(:,4) , day , deltamd(:,5) )
title( 'Daily Mass Deposited' )
xlabel( 'Mission Day' )
ylabel( 'Mass Deposited (g)' )
legend( [ 'View Factor = ' , num2str(F(1)) ] , [ 'View Factor = ' , num2str(F(2)) ] , [ 'View Factor = ' , num2str(F(3)) ] , [ 'View Factor = ' , num2str(F(4)) ] , [ 'View Factor = ' , num2str(F(5)) ] )
end