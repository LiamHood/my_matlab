clear
close all
in2psf = .03613*12^2 ; % inches of water to pounds per square foot
ptare = 1.44 ;
%% Blockage
dp_nozzler = [ 0 .2 .3 .4 .5 ] ; % reading (in H2O)
dp_nozzle = dp_nozzler*in2psf ; % convert to something useful (lb/ft^2)
r = 287.05*5.97994 ; % ft lb/( slug * R )
pa = 14.7*12^2 ; % ambient pressure (lb/ft^2)
ta = 518.47 ; % Rankine
rho = pa/(r*ta) ; % slug/ft^3
rpm = 27 * 58.25 ; %RPM from Hertz

pp = [ .08 .077 .076 .074 .072 ] ; %percent power
powerm = .75*pp ; % convert percent to kW 
power = powerm*737.56 ; % kW to ft lb/s
workrate = (power*(2*pi))/778.16 ; 
d = [ 4 3 1.5 1 .5 ] / 12 ; % diameter in ft
a = (d/2).^2*pi ; % area in ft^2
pfangage = [ 1.68 1.65 1.58 1.57 1.56 ] ; %measuredin H2O
pfan = ( pfangage - ptare )*in2psf ; % useful value lb/ft^2
v = sqrt( ( 2*dp_nozzle ) / rho ) ; % Velocity in ft/s
    %Solve for a reasonable velocity with no blockage by fitting function
    %and extrapolation. I did this because velocity shouldn't be 0
    vfit = polyfit( d(2:5) , v(2:5) , 1 );
    vfun = @(x) vfit(1)*x + vfit(2) ;
efficiency = (a.*v.*(dp_nozzle+pa))./(power) ; % Formula given in lab
mdot = rho.*a.*vfun(d) ; % Mass flow rate
    efit = polyfit( d(2:5) , efficiency(2:5) , 2 );
    efun = @(x) efit(1)*x.^2 + efit(2)*x + efit(3) ;

% Graphs
figure
plot( d , pfan )
title( 'Fan Pressure vs Diameter' )
xlabel( 'Diameter (ft)' )
ylabel( 'Pressure (lb/ft^2)' )
figure
plot( d , dp_nozzle )
title( 'Nozzle Pressure vs Diameter' )
xlabel( 'Diameter (ft)' )
ylabel( 'Pressure (lb/ft^2)' )
figure
plot( d , efficiency )
title( 'Efficiency vs Diameter' )
xlabel( 'Diameter (ft)' )
ylabel( 'Efficiency' )
figure
plot( d , mdot )
title( 'Mass Flow Rate vs Diameter' )
xlabel( 'Diameter (ft)' )
ylabel( 'Mass Flow Rate' )

%% Venturi
pa = linspace( pa , pa , 4 ) ;
pfanvr = [ 1.61 1.53 1.4 1.21 ] ; %measured in H20
dp_nozzvr = [ .25 .45 .7 1.11 ] ; %measured in H20
p1vr = [ 1.58 1.78 2.17 2.56 ] ; %measured in H20
pfanv = ( pfanvr - ptare ) * in2psf ; %dynamic pressure at fan lb/ft^2
dp_nozzv = dp_nozzvr * in2psf ; % dynamic pressure lb/ft^2
p1v = ( p1vr - ptare ) * in2psf ; %dynamic venturi pressure lb/ft^2
ppv = [ 10.7 17.5 27.1 41 ] * .01 ; % percent power
powerv = ppv*.75*737.56 ; % power in ft-lb/s
hz = [ 25 35 45 55 ] ; %Frequency
rpmv = hz*58.25 ; %RPM from Hz
velocityv = sqrt( 2*( dp_nozzv / rho ) ) ; % velocity using dynamic P equation (ft/s)
Q = a(1)*velocityv ; %Volumetric Flow rate
efficiencyv = (Q.*(dp_nozzv+pa))./powerv ; %efficiency using formula. I don't think that pressure is correct since efficiency is over 100 percent

%Graph
figure
plot( rpmv , efficiencyv )
title( 'Efficiency vs rpm' )
xlabel( 'rpm' )
ylabel( 'Efficiency' )
figure
plot( rpmv.^2 , pfanv )
title( 'Fan Pressure vs rpm' )
xlabel( 'rpm' )
ylabel( 'Fan Pressure (lb/ft^2)' )
figure
plot( rpmv.^3 , powerv )
title( 'Power vs rpm' )
xlabel( 'rpm' )
ylabel( 'Power (ft-lb/s)' )
figure
plot( rpmv , velocityv )
title( 'Velocity vs rpm' )
xlabel( 'rpm' )
ylabel( 'Velocity (ft/s)' )

%% Heat
pa = 14.7*12^2 ; %ambient pressure
cp = 187.5387922065*32.174 ; % Cp in pound-force/slug-Rankine
dp_nozzh = .1*in2psf ; % pounds per square foot
pfanhr = 1.68 ; % in H2O
hzh = 27 ; % Hertz
te = 533.07 ; % Rankine
rhoe = (pa+dp_nozzh)/(r*te) ; % Density at exit
I1 = 3.3 ; % Current Amp
I2 = 3.4 ;
V = 120 ; % Voltage (V)
velocityh = sqrt( 2*( dp_nozzh / rhoe ) ) ; % ft/s
mdoth = rhoe*a(1)*velocityh ; %Mass flow rate
    %We don't have percent power for this data. Using polynomials fit to
    %blockage data I made formulas to find what power corresponds to a
    %particular mass flow rate
    dfm = polyfit( mdot(2:5) , d(2:5) , 2 ) ;
    d_from_mdot = @(x) dfm(1)*x.^2 + dfm(2)*x + dfm(3) ;
    pfd = polyfit( power(2:5) , d(2:5) , 2 ) ;
    power_from_d = @(x) pfd(1)*x.^2 + pfd(2)*x + pfd(3) ;
powerh = power_from_d( d_from_mdot( mdoth ) ) ; %Shaft power
aheattran = - powerh + mdoth*( cp*( te - ta ) + ( .5 * velocityh^2 )) ; %Actual heat transfer to air
theattran = .5 * ( I1 + I2 ) * V * .73756 ; %Theoretical heat transfer
efficiencyh = aheattran/theattran ; %Heater efficiency

