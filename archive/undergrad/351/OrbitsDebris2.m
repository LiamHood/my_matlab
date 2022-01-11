%% LEO2LEO

% Sophia Brown, Elliot Rockett, Kenny Carnemolla, Liam Hood
% AERO351 - Orbital Mechanics
% 12/06/18

clc, close, clear all

%% Variable Definiton

mu = 398600; % mu of the Earth
r_earth = 6378; % radius of the Earth

%% Orbiting Debris TLEs to COEs

% Obj1 = GEO Debris (ASIASAT 3S)
% Obj2 = MEO Debris (GIOVE-A)
% Obj3 = LEO Debris #1 (Cosmos 2251 Debris)
% Obj4 = LEO Debris #2 (Jason-1)

Obj11 = '1 25657U 99013A   18330.09763593 -.00000232  00000-0  00000+0 0  9991';
Obj12 = '2 25657   3.5863  75.6443 .0002971 174.2464 356.2930  1.00268800 72145';
        
Obj21 = '1 28922U 05051A   18330.48606153 -.00000043 +00000-0 +00000-0 0  9990';
Obj22 = '2 28922 057.5942 065.8126 .0005092 057.0370 303.0237 01.69405261079977';
    
Obj31 = '1 40796U 93036BVZ 18326.00741122 +.00000423 +00000-0 +88626-4 0  9997';
Obj32 = '2 40796 074.0222 073.7829 .0077119 137.7740 341.2484 14.61190781267082';
    
Obj41 = '1 26997U 01055A   18330.62799917 -.00000074 +00000-0 -49900-4 0  9998';
Obj42 = '2 26997 066.0394 078.0495 .0007059 259.4191 164.6717 12.83995137794334';
    
COE1 = stringTLE_reader(Obj11,Obj12,1) ;
COE2 = stringTLE_reader(Obj21,Obj22,0);
COE3 = stringTLE_reader(Obj31,Obj32,0);
COE4 = stringTLE_reader(Obj41,Obj42,0);
% COE = [epoch,inc,RAAN,ecc,ArgPer,MA,n];

%% GEO Orbit Calculations

n1 = COE1(7)*((2*pi)/(24*60*60)); % Taking mean motion from TLE COE column vector in revs/s
MA1 = COE1(6); % Taking Mean Anomaly from TLE COE column vector 
ecc1 = COE1(4); % Taking eccentricity from TLE COE column vector

T1 = (2*pi)/n1; % Period of GEO orbit in seconds
t1 = MA1/n1; % time since periapsis passage

SMA1 = (mu/(n1)^2)^(1/3); % Semi-major Axis (km)
E1 = eccentric_anomaly(ecc1, MA1); % Calculating Eccentric Anomaly 
TA1 = acos((ecc1 - cos(E1))/(ecc1*cos(E1) - 1)); % Calculating TA in radians
h1 = sqrt(SMA1*mu*(1-ecc1^2)); % Angular Momentum

[r_ppre,v_ppre,r_geopre,v_geopre] = COEs2State(COE1(4), COE1(2), COE1(3), COE1(5), h1, TA1, mu);

state1(1:3) = r_geopre';
state1(4:6) = v_geopre';
t_span1 = [0 T1];
options = odeset('RelTol', 1e-8,'AbsTol', 1e-8);

[t_new1 ,state_new1] = ode45(@OrbitStep, t_span1, state1, options, mu);

% [r_p,v_p,r_geo,v_geo] = COEs2State(ecc, inc, RAAN, ArgPer, h, TA,mu)
% COE = [epoch,inc,RAAN,ecc,ArgPer,MA,n];

% Propogating Orbit Forward
% Start day will be ceil(COE1(1)) = JD2000 = 18331 which equates to
% November 27, 2018
t_prop = ceil(COE1(1)) - COE1(1); % Time to propogate orbit forward in percentage of a day

t_prop = 24*60*60*t_prop; % time to propogate orbit forward in seconds

t_span_move = [0 t_prop];

statemove(1:3) = r_geopre';
statemove(4:6) = v_geopre';
t_span1 = [t_new1(401) t_prop];
options = odeset('RelTol', 1e-8,'AbsTol', 1e-8);
[t_new2 ,state_new1] = ode45(@OrbitStep, t_span1, statemove, options, mu);

% *********** FROM THIS MOMENT ON MISSION IN PROGRESS **************

element = length(state_new1);
% r_geo_start = [state_new1(element,1), state_new1(element,2), state_new1(element,3)];
% v_geo_start = [state_new1(element,4), state_new1(element,5), state_new1(element,6)];

r_geo_start = [state_new1(end,1), state_new1(end,2), state_new1(end,3)]; %these vecs grab all values. element was stopping short
v_geo_start = [state_new1(end,4), state_new1(end,5), state_new1(end,6)];



[h1,inc1,N1,RAAN1,ecc1,ArgPer1,TA1,SMA1,T1,n1] = Orbital_elements(r_geo_start,v_geo_start,mu);

% Start at a true anomaly of 329.41441 degrees and orbit 5 times

Time_since_mission_start = T1*5; % Staying at first orbit for 5 orbits

E = 2*atan((sqrt((1 - ecc1)/(1 + ecc1)))*(tand(TA1/2))); % radians

MA1 = E - ecc1*sin(E);

t = (MA1/(2*pi))*T1; % time since periapsis passage

t_add = T1 - t; % time to get back to periapsis
t_add1 = T1*0.5; % time to get to apoapsis for inc change

t_span1 = [0 t_add];
[t_new ,state_new1] = ode45(@OrbitStep, t_span1, state1, options, mu);
%Geo plot
figure
view(43,24);
[x,y,z] = sphere;
surf(6378*x,6378*y,6378*z);
axis([-5e+4 5e+4 -5e+4 5e+4 -5e+4 5e+4]);
xlabel('x, km');
ylabel('y, km');
zlabel('z, km');
%Animated line plot either use plot3 or code below
% curve = animatedline('LineWidth',2);
% for i = 1:length(state_new1)
%     addpoints(curve,state_new1(i,1),state_new1(i,2),state_new1(i,3));
%      drawnow
%     pause(0.001);
% end
 hold on 
plot3(state_new1(:,1), state_new1(:,2), state_new1(:,3));
plot3(state_new1(1,1), state_new1(1,2), state_new1(1,3), '*') %plot start



Time_since_mission_start = Time_since_mission_start + t_add + t_add1;
% This is the time for 5 orbits, getting back to periapsis and then
% continuing on to apoapsis

% Inc change and start of transfer orbit

n2 = COE2(7)*((2*pi)/(24*60*60)); % Taking mean motion from TLE COE column vector in revs/s
MA2 = COE2(6); % Taking Mean Anomaly from TLE COE column vector 
ecc2 = COE2(4); % Taking eccentricity from TLE COE column vector

T2 = (2*pi)/n2; % Period of MEO orbit in seconds
t2 = MA2/n2; % time since periapsis passage

SMA2 = (mu/(n2)^2)^(1/3); % Semi-major Axis (km)
E2 = eccentric_anomaly(ecc2, MA2); % Calculating Eccentric Anomaly 
TA2 = acos((ecc2 - cos(E2))/(ecc2*cos(E2) - 1)); % Calculating TA in radians
h2 = sqrt(SMA2*mu*(1-ecc2^2)); % Angular Momentum
[r_ptran,v_ptran,r_geotran,v_geotran] = COEs2State(COE2(4), COE2(2), COE2(3), COE2(5), h2, TA2, mu);

t_span1 = [t_new1(401) T2];
state_new2 = [r_ptran, v_ptran];
[t_new1 ,state_new3] = ode45(@OrbitStep, t_span1, state_new2, options, mu);
plot3(state_new3(:,1), state_new3(:,2), state_new3(:,3), 'r')

[h2,inc2,N2,RAAN2,ecc2,ArgPer2,TA2,SMA2,T2,n2] = Orbital_elements(r_geotran,v_geotran,mu);

orbit2_periapsis = SMA2*(1 - ecc2);
orbit1_apoapsis = SMA1*(1 + ecc1);

ecctran = (orbit1_apoapsis - orbit2_periapsis)/(orbit2_periapsis + orbit1_apoapsis);
htran = sqrt(orbit1_apoapsis*mu*(1 + ecctran));

va_orbit1 = (mu/h1)*(1 + ecc1);
va_tran = (mu/h2)*(1 - ecctran);
vp_tran = (mu/h2)*(1 + ecctran);
SMAtran = (orbit1_apoapsis + orbit2_periapsis)/2;
deltainc = norm(COE1(2) - COE2(2));

deltav1 = sqrt(va_orbit1^2 + va_tran^2 - 2*va_orbit1*va_tran*cos(deltainc)); % changing inc and getting onto transfer orbit

                                % Orbital elements and state of transfer
                                %(ecc, inc, RAAN, ArgPer, h, TA, mu)
                                incTrans1 = COE2(2) ;       
                                [~,~,r_Trans1,v_Trans1] = COEs2State(ecctran, incTrans1, RAAN1, ArgPer1 , htran, 0, mu) ;
        
       
Ttran = ((2*pi)/sqrt(mu))*SMAtran^(3/2);
t_add2 = Ttran*0.5;

%plot trans orbit of inc change
state_trans = [r_ptran,v_ptran];
t_span = [t_new1(265) Ttran];
[t_new2 ,state_new1] = ode45(@OrbitStep, t_span, state_trans, options, mu);
hold on
plot3(state_new1(:,1), state_new1(:,2), state_new1(:,3) ,'r')

Time_since_mission_start = Time_since_mission_start + t_add2;
% Includes the time for 5 orbits, getting back to periapsis,
% continuing on to apoapsis and traveling to perigee of transfer orbit


    % ******* NEED TO TALK ABOUT APSE LINE CHANGE IF NEEDED *******
    % Going to talk to DR. A on Monday to answer questions
    
    
%% MEO Orbit Calculations

vp_orbit2 = (mu/h2)*(1 + ecc2); % velocity at perigee of orbit 2
deltav2 = vp_tran - vp_orbit2 ; % delta v burn to get off transfer orbit and onto MEO orbit
deltav_tot = deltav1 + deltav2 ; % running total for our delta v (KEEP ADDING TO THIS)

% finding location of MEO object at this time

t_meo = ceil(COE1(1)) + (Time_since_mission_start/(3600*24)) ; % finding the JD2000 time (days) at this point in mission in order to propogate MEO satellite position
% Propogating Orbit Forward
t_prop2 = t_meo - COE2(1); % Time to propogate orbit forward in percentage of a day
t_prop2 = 24*60*60*t_prop2; % time to propogate orbit forward in seconds
t_span_move2 = [0 t_prop2];

[r_p2,v_p2,r_geo2,v_geo2] = COEs2State(COE2(4), COE2(2), COE2(3), COE2(5), h2, TA2, mu); % getting state vector of MEO object from TLE

statemove2(1:3) = r_geo2';
statemove2(4:6) = v_geo2';
options = odeset('RelTol', 1e-8,'AbsTol', 1e-8);
[t_new4 ,state_new2] = ode45(@OrbitStep, t_span_move2, statemove2, options, mu); % getting new state vector of MEO object at this time in mission 


element2 = length(state_new2);
% r_meo = [state_new2(element2,1), state_new2(element2,2), state_new2(element2,3)]; % position vector of meo object
% v_meo = [state_new2(element2,4), state_new2(element2,5), state_new2(element2,6)]; % velocity vector of meo object
r_meo = [state_new2(end,1), state_new2(end,2), state_new2(end,3)]; % position vector of meo object
v_meo = [state_new2(end,4), state_new2(end,5), state_new2(end,6)]; % velocity vector of meo object
t_span_move2 = [t_new2(129) t_prop2];
statemove2 = [r_meo, v_meo];
[t_new2 ,state_new2] = ode45(@OrbitStep, t_span_move2, statemove2, options, mu); % getting new state vector of MEO object at this time in mission 
 %Meo plot
hold on
plot3(state_new2(:,1), state_new2(:,2), state_new2(:,3));

 % getting coes of MEO orbit right now
[hmeo,incmeo,Nmeo,RAANmeo,eccmeo,ArgPermeo,TAmeo,SMAmeo,Tmeo,nmeo] = Orbital_elements(r_meo,v_meo,mu) ; 

% Using these COEs, the MEO object is at a true anomaly of 265.1025
% degrees, and our spacecraft is at perigee (true anomaly of 0). We need to
% do a phasing maneuver so that they intercept

% outer transfer phasing maneuver since the chaser(spacecraft) is ahead of
% the target(satellite)
TAmeo = deg2rad(360 - TAmeo) ; % converting true anomaly from degrees to radians
Emeo = 2*atan(sqrt((1-eccmeo)/(1+eccmeo))*tan(TAmeo/2)); % eccentric anomaly for MEO orbit
Memeo = Emeo - eccmeo*sin(Emeo) ; % mean anomaly for MEO orbit

t_meo_per = (Memeo/nmeo) ; % the time it will take the MEO satellite to get from current true anomaly to perigee
tphase_meo = Tmeo + t_meo_per ; % time of the phase orbit so both objects reach perigee at same time

aphase_meo = ((tphase_meo*sqrt(mu))/(2*pi))^(2/3) ; % semimajor axis of the phase orbit
ra_phasemeo = (2*aphase_meo) - orbit2_periapsis ; % apoapse of the phase orbit
eccphase_meo = (ra_phasemeo - orbit2_periapsis)/(orbit2_periapsis + ra_phasemeo) ; % eccentricity of the phase orbit
hphase_meo = sqrt(ra_phasemeo*mu*(1+eccphase_meo)) ; % specific angular momentum of the phase orbit

                % COES and state of phasing
                [~,~,r_Phase1,v_Phase1] = COEs2State(eccphase_meo, incmeo, RAANmeo, ArgPermeo, hphase_meo, 0, mu) ;
                                        
dv_meo1 = (hphase_meo/orbit2_periapsis) - vp_orbit2 ; % change in velocity to get onto phase orbit
dv_meo2 = vp_orbit2 - (hphase_meo/orbit2_periapsis) ; % change in velocity to get back onto MEO orbit

dv_phase = abs(dv_meo2) + abs(dv_meo1) ; % delta v for phase maneuver

deltav_tot = deltav_tot + dv_phase ; % new delta v total

% getting state vector of spacecraft at perigee in MEO orbit
t_span_meo = [t_new2(129) t_meo_per];

statemeo2(1:3) = r_meo';
statemeo2(4:6) = v_meo';
options = odeset('RelTol', 1e-8,'AbsTol', 1e-8);
[t_new5 ,state_meo2] = ode45(@OrbitStep, t_span_meo, statemeo2, options, mu); % getting new state vector of MEO object at this time in mission 
%Meo plot
hold on
plot3(state_meo2(:,1), state_meo2(:,2), state_meo2(:,3))
%Animated line plot either use plot3 or code below
% curve = animatedline('LineWidth',1,'Color','r');
% for i = 1:length(state_meo2)
%     addpoints(curve,state_meo2(i,1),state_meo2(i,2),state_meo2(i,3));
%      drawnow
%     pause(0.001);
% end
r_meop = [state_meo2(1), state_meo2(2), state_meo2(3)]; % position vector of meo object
v_meop = [state_meo2(4), state_meo2(5), state_meo2(6)]; % velocity vector of meo object

Time_since_mission_start = Time_since_mission_start + tphase_meo ; % new time since mission start

% at this point, the MEO object and the spacecraft are both at perigee for
% the MEO orbit. Now they will orbit 5 times

Time_since_mission_start = Time_since_mission_start + (Tmeo*5) ;

%% Now the transfer to the first LEO object : Cosmos 2251 Debris

% bi-elliptic hohmann transfer - find state vector for object 3 at this time.
t_leo1 = ceil(COE1(1)) + (Time_since_mission_start/(3600*24)) ; % finding the JD2000 time (days) at this point in mission in order to propogate MEO satellite position

% Propogating Orbit Forward

t_prop3 = t_leo1 - COE3(1); % Time to propogate orbit forward in percentage of a day

t_prop3 = 24*60*60*t_prop3; % time to propogate orbit forward in seconds

t_span_move3 = [0 t_prop3];

n3 = COE3(7)*((2*pi)/(24*60*60)); % Taking mean motion from TLE COE column vector in revs/s
MA3 = COE3(6); % Taking Mean Anomaly from TLE COE column vector 
ecc3 = COE3(4); % Taking eccentricity from TLE COE column vector

T3 = (2*pi)/n3; % Period of LEO orbit in seconds
t3 = MA1/n3; % time since periapsis passage

SMA3 = (mu/(n3)^2)^(1/3); % Semi-major Axis (km)
E3 = eccentric_anomaly(ecc3, MA3); % Calculating Eccentric Anomaly 
TA3 = acos((ecc3 - cos(E3))/(ecc3*cos(E3) - 1)); % Calculating TA in radians
h3 = sqrt(SMA3*mu*(1-ecc3^2)); % Angular Momentum

[r_p3,v_p3,r_geo3,v_geo3] = COEs2State(COE3(4), COE3(2), COE3(3), COE3(5), h3, TA3, mu); % getting state vector of LEO object from TLE

statemove3(1:3) = r_geo3';
statemove3(4:6) = v_geo3';
options = odeset('RelTol', 1e-8,'AbsTol', 1e-8);
[t_new6 ,state_new3] = ode45(@OrbitStep, t_span_move3, statemove3, options, mu); % getting new state vector of LEO object at this time in mission 


element3 = length(state_new3);
r_leo1 = [state_new3(element3,1), state_new3(element3,2), state_new3(element3,3)]; % position vector of leo object
v_leo1 = [state_new3(element3,4), state_new3(element3,5), state_new3(element3,6)]; % velocity vector of leo object

% getting coes of LEO orbit right now
[h3,inc3,N3,RAAN3,ecc3,ArgPer3,TA3,SMA3,T3,n3] = Orbital_elements(r_leo1,v_leo1,mu) ; 

rp_leo1 = ((h3^2)/mu)*(1/(1+ecc3)) ; % rp of leo1 orbit
ra_leo1 = ((h3^2)/mu)*(1/(1-ecc3)) ; % ra of leo1 orbit
state_leo1 = [r_leo1,v_leo1];
t_span = [0 T3];
[t_new7 ,state_leo1] = ode45(@OrbitStep, t_span, state_leo1, options, mu);
 hold on
 plot3(state_leo1(:,1), state_leo1(:,2), state_leo1(:,3))
 %Animated line plot either use plot3 or code below
%curve = animatedline('LineWidth',2,'Color','b');
% for i = 1:length(state_leo1)
%     addpoints(curve,state_leo1(i,1),state_leo1(i,2),state_leo1(i,3));
%      drawnow
%     pause(0.001);
% end

ra_meo = ((hmeo^2)/mu)*(1/(1-eccmeo)) ; % ra of meo orbit
ra_ellipse = 13*ra_leo1 ;

e_trans1 = (ra_ellipse - orbit2_periapsis)/(orbit2_periapsis + ra_ellipse) ; % eccentricity of the first transfer ellipse
a_trans1 = (orbit2_periapsis + ra_ellipse)/2 ; % semi major axis of first transfer ellipse
T_trans1 = ((2*pi)/(sqrt(mu)))*(a_trans1^(3/2)) ; % period of first transfer ellipse
h_trans1 = sqrt(ra_ellipse*mu*(1-e_trans1)) ; % specific angular momentum of first transfer ellipse


                % COEs and state of first ellipse
                TA_trans1 = 0 ;
                [~,~,r_ellipse1,v_ellipse1] = COEs2State(e_trans1, inc2, RAAN2, ArgPer2, h_trans1, TA_trans1, mu) ;
                state_ellipse1o = [ r_ellipse1 ; v_ellipse1 ] ;
                t_ellipse1 = [ 0 , T_trans1/2 ];
                [t_ellipse1 , state_ellipse1] = ode45(@OrbitStep, t_ellipse1, state_ellipse1o, options, mu);                    
                                    
e_trans2 = (ra_ellipse - rp_leo1)/(ra_ellipse + rp_leo1) ; % eccentricity of the 2nd transfer ellipse
a_trans2 = (rp_leo1 + ra_ellipse)/2 ; % semi major axis of 2nd transfer ellipse
T_trans2 = ((2*pi)/(sqrt(mu)))*(a_trans2^(3/2)) ; % period of 2nd transfer ellipse
h_trans2 = sqrt(ra_ellipse*mu*(1-e_trans2)) ; % specific angular momentum of 2nd transfer ellipse

                % COEs and state of first ellipse
                TA_trans2 = 180 ;
                [~,~,r_ellipse2,v_ellipse2] = COEs2State(e_trans2, inc3, RAAN2, ArgPer2, h_trans2, TA_trans2, mu) ;
                state_ellipse2o = [ r_ellipse2 ; v_ellipse2 ] ;
                t_ellipse2 = [ 0 , T_trans2/2 ];
                [t_ellipse1 , state_ellipse2] = ode45(@OrbitStep, t_ellipse2, state_ellipse2o, options, mu);                  

                
                figure
                hold on 
                plot3( state_ellipse1(:,1) , state_ellipse1(:,2) , state_ellipse1(:,3) )
                plot3( state_ellipse2(:,1) , state_ellipse2(:,2) , state_ellipse2(:,3) )
                hold off
                
%             v1_trans1 = h_trans1/orbit2_periapsis ; % perigee of first ell
%             v1_trans2 = h_trans2/rp_trans1 ; % perigee of second ellipse
%             v2_trans2 = h3/rp_leo1 ; % perigee speed of first leo
%             deltav_trans = abs(vp_orbit2 - v1_trans1) + abs(v1_trans1 - v1_trans2) + abs(v1_trans2 - v2_trans2) ; 

vp_trans1 = h_trans1/orbit2_periapsis ; % perigee of first ellipse
va_trans1 = h_trans1/ra_ellipse ;
vp_trans2 = h_trans2/rp_leo1 ; % perigee of second ellipse
va_trans2 = h_trans2/ra_ellipse ;
vp3 = h3/rp_leo1 ; % perigee speed of first leo
deltav_incM2L = sqrt( va_trans1 + va_trans2 - 2*va_trans1*va_trans2*cos( abs( COE2(2) - COE3(2) ) ) ) ;
deltav_trans = abs( vp_orbit2 - vp_trans1) + deltav_incM2L + abs( vp_trans2 - vp3 ) ; 

deltav_tot = deltav_tot + deltav_trans; % new deltav_total

t_trans1 = 0.5*T_trans1 ; % time for transfer 1
t_trans2 = 0.5*T_trans2 ; % time for transfer 2

Time_since_mission_start = Time_since_mission_start + t_trans1 + t_trans2 ; 

%% Now, spacecraft is at perigee of COSMOS 2251 debris. Need to use 
% Time_since_mission_start to figure out where Cosmos debris is located in
% the orbit right now and then do a phasing maneuver so that the spacecraft
% and the debris intersect.

%% LEO2LEO
    % find position after transfer to its orbit
    % phase to match positions
    % lamberts from first to where ever second LEO object is after transfer
    % use smallest delta v lamberts 
    
% first LEO Object
        n3 = COE3(7)*((2*pi)/(24*60*60)); % Taking mean motion from TLE COE column vector in revs/s
        MA3 = COE3(6); % Taking Mean Anomaly from TLE COE column vector 
        ecc3 = COE3(4); % Taking eccentricity from TLE COE column vector

        T3 = (2*pi)/n3; % Period of GEO orbit in seconds
        t3 = MA3/n3; % time since periapsis passage
        t_afterM2L = Time_since_mission_start + 18331*60*60*24 - COE3(1)*60*60*24 ; % time from epoch to arrival on first LEO orbit
        t_fromPer_beg = (t_afterM2L + t3) - T3*floor((t_afterM2L + t3)/T3) ; % time LEO object is from perigee at arrival from MEO
        MA3 = t_fromPer_beg * n3 ;
        SMA3 = (mu/(n3)^2)^(1/3); % Semi-major Axis (km)
        E3 = eccentric_anomaly(ecc3, MA3); % Calculating Eccentric Anomaly 
        TA3 = 2*pi + 2*atan( tan(E3/2) / sqrt((1-ecc3)/(1+ecc3)) ); % Calculating TA in radians
        h3 = sqrt(SMA3*mu*(1-ecc3^2)); % Angular Momentum
        rp3 = (h3^2/mu)*(1/(1+ecc3)) ;
        vp3 = mu/h3*(1+ecc3) ;

    % Phase to first LEO sat
        % Outer because target is trailing
        T_phase3 = T3 + t_fromPer_beg ; % period of phasing orbit to reach object 3
        aPhase3 = ( ( T_phase3 * sqrt(mu) )/(2*pi) )^(2/3) ; % semi-major of phasing orbit
        ra3phase = 2*aPhase3 - rp3 ;
        ecc3phase = (-rp3+ra3phase)/(rp3+ra3phase) ; % perigee of first LEO is apogee of phasing orbit
        h3phase = sqrt( rp3*mu*(1+ecc3phase) ) ;
        vp3phase = mu/h3phase*(1+ecc3phase) ;
        
                            % COEs and state of first ellipse
                            TA3phase = 0 ;
                            [~,~,r_phase2,v_phase2] = COEs2State(ecc3phase, COE3(2), COE3(3), COE3(5), h3phase, TA3phase, mu) ;
                                                            
        dv3phase = (vp3 - vp3phase)*2 ; % dv from phasing
        deltav_tot = abs(dv3phase) + deltav_tot ; % update total
        Time_since_mission_start = Time_since_mission_start + T_phase3 + T3*5 ;
        [~,~,r_3,v_3] = COEs2State(COE3(4), COE3(2), COE3(3), COE3(5) , h3 , 0 , mu);
state_trans = [r_3,v_3];
t_span = [0 T_phase3];
[t_new8 ,state_new1] = ode45(@OrbitStep, t_span, state_trans, options, mu);
hold on
plot3(state_new1(:,1), state_new1(:,2), state_new1(:,3))
    % second LEO
        n4 = COE4(7)*((2*pi)/(24*60*60)); % Taking mean motion from TLE COE column vector in revs/s
        MA4 = COE4(6); % Taking Mean Anomaly from TLE COE column vector 
        ecc4 = COE4(4); % Taking eccentricity from TLE COE column vector
        T4 = (2*pi)/n4; % Period of GEO orbit in seconds
        
    %transfer from first to second LEO object     
         t_L2Lall = linspace( 1e4 , 3e4 , 5e2 ); % time to transfer from LEO to LEO
        for ii = 1:length(t_L2Lall)
            progress = ii/length(t_L2Lall)*100
            t_L2L = t_L2Lall(ii) ;
            t4 = MA4/n4; % time since periapsis passage
            t_afterL2L = Time_since_mission_start + t_L2L + 18331*60*60*24 - COE4(1)*60*60*24 ; % time from epoch to arrival on first LEO orbit
            t_fromPer_beg4 = (t_afterL2L + t4) - T4*floor((t_afterL2L + t3)/T4) ; % time LEO object is from perigee at arrival from MEO
            MA4 = t_fromPer_beg4 * n4 ;
            a4 = (mu/(n4)^2)^(1/3); % Semi-major Axis (km)
            E4 = eccentric_anomaly(ecc4, MA4); % Calculating Eccentric Anomaly 
            TA4 = 2*pi + 2*atan( tan(E4/2) / sqrt((1-ecc4)/(1+ecc4)) ); % Calculating TA in radians
            h4 = sqrt(a4*mu*(1-ecc4^2)); % Angular Momentum
       % Use Lamberts to transfer
            [~,~,r_4,v_4] = COEs2State(COE4(4), COE4(2), COE4(3), COE4(5) , h4 , TA4 , mu); % find state of 4th object after transfer time
            [ v1 , v2  ] = Lamberts2( r_3 , r_4 , t_L2L , mu , 1e-5 , 1 ) ;
            dvBeg = norm( v1 - v_3 ) ;
            dvEnd = norm( v2 - v_4 ) ;
            dv(ii) = dvBeg+dvEnd ;
        end
        
          % plot trans orbit
state_trans = [r_4,v_4];
t_span = [0 t_L2Lall];
[t_new9 ,state_new1] = ode45(@OrbitStep, t_span, state_trans, options, mu);
hold on
plot3(state_new1(:,1), state_new1(:,2), state_new1(:,3))
        % Optimization
        figure
            plot( t_L2Lall , dv )
            ylabel( 'Delta v (km/s)' )
            xlabel( 'Transfer Time (s)' )
            [ dvOpt , i ] = min( dv ) ;
        deltav_tot = deltav_tot + dvOpt ;
        Time_since_mission_start = Time_since_mission_start + t_L2Lall(i) ;
        MissionTime = Time_since_mission_start/(60*60*24) ;
        
%% Functions

function [COE] = stringTLE_reader(TLE1, TLE2, def)
%This function takes a TLE inputed as 2 strings and outputs the COEs. The
%def defines whether there is a character space between the mean motion
%and obrit # which will determine how the functions grabs the COE. 1 = yes
%0 = no

A = strsplit(TLE1,' ');
B = strsplit(TLE2,' ');

epoch = str2double(A{4}); % Epoch date in Julian format
inc = deg2rad(str2double(B{3})); % Inclination in radians
RAAN = deg2rad(str2double(B{4})); % RAAN in radians
ecc = str2double(B{5}); % Eccentricity
ArgPer = deg2rad(str2double(B{6})); % Argument of Perigee in radians
MA = deg2rad(str2double(B{7})); % Mean Anomaly in radians

if def == 1 % Mean Motion 
    
    n = str2double(B{8});
    
else
    
    a = B{8};
    n = str2double(a(1:11));
    
    
end

COE = [epoch,inc,RAAN,ecc,ArgPer,MA,n];

    
end

function [E] = eccentric_anomaly(ecc, MeA)
% This function uses Kepler's equation to find the eccentric anomaly given 
% the eccentricity and mean anomaly

ratio_limit = 1e-8;
ratio = 1;

if MeA < pi
    
    E = MeA + ecc/2;
    
else
    
    E = MeA - ecc/2;
    
end

while abs(ratio) > ratio_limit
    
f = E - ecc*sin(E) - MeA;
fprime = 1 - ecc*cos(E);

ratio = f/fprime;

E = E - ratio;
   
end

end

function [r_p,v_p,r_geo,v_geo] = COEs2State(ecc, inc, RAAN, ArgPer, h, TA, mu)
% This function takes 6 classical orbital elements and computes the state
% vector, ie the position (r) and velocity (v) vectors.

% r_p = r vector in perifocal frame
% v_p = v vector in perifocal frame

% r_geo = r vector in geocentric frame
% v_geo = v vector in geocentric frame

r_p(1,1) = (h^2/mu)*(1/(1+(ecc*cos(TA))))*(cos(TA));
r_p(2,1) = (h^2/mu)*(1/(1+(ecc*cos(TA))))*(sin(TA));
r_p(3,1) = 0;

v_p(1,1) = (mu/h)*-sin(TA);
v_p(2,1) = (mu/h)*(ecc + cos(TA));
v_p(3,1) = 0;

Qx = [-sin(RAAN)*cos(inc)*sin(ArgPer) + cos(RAAN)*cos(ArgPer), -sin(RAAN)*cos(inc)*cos(ArgPer) - cos(RAAN)*sin(ArgPer), sin(RAAN)*sin(inc);
      cos(RAAN)*cos(inc)*sin(ArgPer) + sin(RAAN)*cos(ArgPer), cos(RAAN)*cos(inc)*cos(ArgPer) - sin(RAAN)*sin(ArgPer), -cos(RAAN)*sin(inc);
      sin(inc)*sin(ArgPer), sin(inc)*cos(ArgPer), cos(inc)];
  
r_geo = Qx*r_p;
v_geo = Qx*v_p;

end

function [h,inc,N,RAAN,ecc,ArgPer,TA,SMA,T,n] = Orbital_elements(r,v,mu)
% Given a position and velocity vector, this function calculates the
% following orbital elements

% h - angular momentum
% inc - inclination
% N - node line
% RAAN - Right Ascension of the Ascending Node
% ecc - eccentricity
% ArgPer - argument of perigee
% TA = true anomaly
% SMA = semi-major axis
% T = Period
% n = mean motion

K = [0,0,1];

rmag = sqrt(dot(r,r));
v_radial = dot(r,v)/rmag;

hvec = cross(r,v);
h = sqrt(dot(hvec,hvec));

inc = acosd(hvec(3)/h);

Nvec = cross(K,hvec);
N = sqrt(dot(Nvec,Nvec));

if Nvec(2) < 0
    
    RAAN = 360 - acosd(Nvec(1)/N);
    
else
    
     RAAN = acosd(Nvec(1)/N);
     
end

eccvec = (1/mu)*(cross(v,hvec) - (mu*(r/rmag)));
ecc = sqrt(dot(eccvec,eccvec));

if eccvec(3) < 0
    
    ArgPer = 360 - acosd((dot(Nvec,eccvec))/(N*ecc));
    
else
    
    ArgPer = acosd((dot(Nvec,eccvec))/(N*ecc));
    
end

if v_radial < 0
    
    TA = 360 - acosd((1/ecc)*((h^2/(mu*rmag)) - 1));
    
else
    
    TA = acosd((1/ecc)*((h^2/(mu*rmag)) - 1));
    
end

SMA = (h^2/mu)*(1/(1-ecc^2));

T = ((2*pi)/sqrt(mu))*SMA^(3/2);

n = (2*pi)/T;

end

function dstatedt = OrbitStep(t,state,muearth)
        
    dx = state(4);
    dy = state(5);
    dz = state(6);
        
    r = norm([state(1) state(2) state(3)]);
        
    ddx = -muearth*state(1)/r^3;
    ddy = -muearth*state(2)/r^3;
    ddz = -muearth*state(3)/r^3;
        
    dstatedt = [dx;dy;dz;ddx;ddy;ddz];  
        
end

function [ v1_long , v2_long , v1_short , v2_short ] = Lamberts( r1 , r2 , dt , mu , tol , pro )
% pro is 1 or 0 for prograde or retrograde respectively
    
    r1mag = norm( r1 ) ;
    r2mag = norm( r2 ) ;
    rcross = cross( r1 , r2 ) ;
    
    % Find delta theta
        if pro == 1 
            if rcross(3) >= 0
                dtheta = acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            else
                dtheta = 2*pi - acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            end
        else
            if rcross(3) < 0 
                dtheta = acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            else
                dtheta = 2*pi - acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            end
        end  
    tm = [ -1 , 1 ] ;
    for ii = 1:2
    dtheta = asin( tm(ii)*(1-cos(dtheta)^2) ) ;
    A = sin( dtheta )*sqrt( r1mag*r2mag/(1-cos(dtheta)) ) ;
        z = 0 ;
        C = 1/2 ;
        S = 1/6 ;
        zup = 4*pi^2 ;
        zlow = -4*pi^2 ;
        y = r1mag + r2mag + (A*(z*S-1))/sqrt(C) ;
        chi = sqrt(y/C) ;
        dtloop = (chi^3*S)/sqrt(mu) + (A*sqrt(y))/sqrt(mu) ;
        while abs( dtloop - dt ) > tol
                if dtloop <= dt
                    zlow = z ;
                else 
                    zup = z ;
                end
             z = ( zup + zlow ) / 2 ;
            [ S , C ] = Stumpf( z ) ;
            y = r1mag + r2mag + (A*(z*S-1))/sqrt(C) ;
            chi = sqrt(y/C) ;
            dtloop = (chi^3*S)/sqrt(mu) + (A*sqrt(y))/sqrt(mu) ;
        end
        f = 1 - y/r1mag ;
        g = A*sqrt(y/mu) ;
        gdot = 1 - y/r2mag ;
        
        v1(:,ii) = ( 1/g )*( r2 - f*r1 ) ;
        v2(:,ii) = ( 1/g )*( gdot*r2 - r1 ) ;
    end
    v1_long = v1(:,1) ;
    v2_long = v2(:,1) ;
    v1_short = v1(:,2) ;
    v2_short = v2(:,2) ;
end

function [ v1 , v2 ] = Lamberts2( r1 , r2 , dt , mu , tol , pro )
% pro is 1 or 0 for prograde or retrograde respectively
    
    r1mag = norm( r1 ) ;
    r2mag = norm( r2 ) ;
    rcross = cross( r1 , r2 ) ;
    
    % Find delta theta
        if pro == 1 
            if rcross(3) >= 0
                dtheta = acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            else
                dtheta = 2*pi - acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            end
        else
            if rcross(3) < 0 
                dtheta = acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            else
                dtheta = 2*pi - acos( dot(r1,r2)/(r1mag*r2mag) ) ;
            end
        end  
    
    A = sin( dtheta )*sqrt( r1mag*r2mag/(1-cos(dtheta)) ) ;
        z = 0 ;
        C = 1/2 ;
        S = 1/6 ;
        zup = 4*pi^2 ;
        zlow = -4*pi^2 ;
        y = r1mag + r2mag + (A*(z*S-1))/sqrt(C) ;
        chi = sqrt(y/C) ;
        dtloop = (chi^3*S)/sqrt(mu) + (A*sqrt(y))/sqrt(mu) ;
        while abs( dtloop - dt ) > tol
                if dtloop <= dt
                    zlow = z ;
                else 
                    zup = z ;
                end
             z = ( zup + zlow ) / 2 ;
            [ S , C ] = Stumpf( z ) ;
            y = r1mag + r2mag + (A*(z*S-1))/sqrt(C) ;
            chi = sqrt(y/C) ;
            dtloop = (chi^3*S)/sqrt(mu) + (A*sqrt(y))/sqrt(mu) ;
        end
        f = 1 - y/r1mag ;
        g = A*sqrt(y/mu) ;
        gdot = 1 - y/r2mag ;
        
        v1 = ( 1/g )*( r2 - f*r1 ) ;
        v2 = ( 1/g )*( gdot*r2 - r1 ) ;
end

function [ S , C ] = Stumpf( z ) 
% Stumpff Functions
sl = 11 ;
sc = zeros(1,sl) ;
cc = zeros(1,sl) ;
S = 0 ;
C = 0 ;

        for kk = 1:sl
            sc(kk) = (-1)^(kk-1) / factorial(2*(kk-1)+3) ;
        end
        for kk = 1:sl
            cc(kk) = (-1)^(kk-1) / factorial(2*(kk-1)+2) ;
        end
        
        for jj = 1:sl
            Scomp(jj) = sc(jj)*z^(jj-1) ;
        end
        S = sum(Scomp) ;
        for jj = 1:sl
            Ccomp(jj) = cc(jj)*z^(jj-1) ;
        end
        C = sum(Ccomp) ;
end