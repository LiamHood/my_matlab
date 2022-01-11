%% Given Values
    r_e = 6378 ; %Radius of Earth in km
    r_v = 6052 ; %Radius of Venus in km
    a_eo = 149.6 * 10^6 ; %Semimajor axis of Earth's orbit in km
    a_vo = 108.2 * 10^6 ; %Semimajor axis of Venus's orbit in km
    mu_e = 398600 ; %Mu of Earth
    mu_v = 324900 ; %Mu of Venus
    mu_s = 1.327 * 10^11 ; %Mu of Sun
    g = 9.81 ; %Earth's Gravity
    a_ep = 300 + r_e ; %semimajor axis of parking orbit of Earth
    a_vp = 1000 + r_v ; %semimajor axis of parking orbit of Venus
    
%% Part A
  disp( 'Part a' )
    
    %Heliocentric
    E_t = - mu_s / ( a_eo + a_vo ) ; %Specific energy of heliocentric transfer orbit
    v_et = ( 2 * (( mu_s / a_eo ) + E_t )) ^.5 ; %Velocity of hto at Earth
    v_vt = ( 2 * (( mu_s / a_vo ) + E_t )) ^.5 ; %Velocity of hto at Venus
    v_e = ( 2 * (( mu_s / a_eo ) - ( mu_s / ( 2 * a_eo )))) ^.5 ; %Velocity of Earth
    v_v = ( 2 * (( mu_s / a_vo ) - ( mu_s / ( 2 * a_vo )))) ^.5 ; %Velocity of Earth
    v_ee = v_et - v_e ; %Velocity of HOOPLAH escaping Earth's influence
    v_ve = v_vt - v_v ; %Velocity of HOOPLAH escaping Venus's influence

    %Within the influence of a planet
    E_ee = ( v_ee ^ 2 ) / 2 ; %Specific energy of orbit to leave Earth's influence
    E_ve = ( v_ve ^ 2 ) / 2 ; %Specific energy of orbit to enter Venus's influence
    E_ep = - mu_e / ( 2 * a_ep ) ; %Specific energy of parking orbit around Earth
    E_vp = - mu_v / ( 2 * a_vp ) ; %Specific energy of parking orbit around Venus
    v_ep = ( 2 * ( mu_e / a_ep + E_ep )) ^ .5 ; %velocity of parking orbit around Earth
    v_vp = ( 2 * ( mu_v / a_vp + E_vp )) ^ .5 ; %velocity of parking orbit around Venus
    v_epe = ( 2 * ( mu_e / a_ep + E_ee )) ^ .5 ; %velocity at beginning of hyberbolic orbit to leave Earth parking orbit
    v_vpe = ( 2 * ( mu_v / a_vp + E_ve )) ^ .5 ; %velocity at end of hyberbolic orbit to arrive at Venus parking orbit

    %delta V calculations
    dve = v_epe - v_ep ; %delta v to leave Earth
    dvv = v_vpe - v_vp ; %delta v to arrive at Venus
    dv = dve + dvv ; %total delta v 
    disp( [ 'HOOPLAH requires a total delta V of ' , num2str( dv ) , ' to go from Earth to Venus' ] )
    
%% Part B
disp( 'Part b' )
    T = 900 * 10^3 ; %Thrust in neutons
    mdot = 300 ; %mass flow rate in kg/s
    m = 10000 ; %mass in kg
    x = dv*mdot;
       
    disp( 'I did not have the equations I needed for this part' )