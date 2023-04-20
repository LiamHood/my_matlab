function [ t , states] = encke_tether( tspan , sc_state0, tether_state0, tether_param, mu , tol )
% Uses Classical Orbital Elements a (semimajor axis), e (eccentricity), 
% i (inclination), RAAN (right ascension of ascending node), aop (argument
% of periapsis), ta (true anomaly)
%
% Considers electrodynamic force and gravity gradient
% will at some point maybe: atmospheric drag, spherical harmonics,
% third body of sun and moon, srp
% tether is treated as dumbbell model, rigid with lumped masses at the end


%     opts = odeset('RelTol', tol, 'AbsTol', tol ) ;
%     [ t , states ] = ode45(@gauss_variations, tspan , ...
%         [ sc_state0 ; tether_state0], opts, tether_param, mu) ;
% 
% 
%     function dstate = gauss_variations(t, states, tether_param, mu)
    dt = 1;
    n = (tspan(2)-tspan(1))/dt;
    
    states = [ sc_state0 ; tether_state0];
    %set size of results
    r = zeros(3,n);
    v = zeros(3,n);
    t = zeros(1,n);
    phi = zeros(1,n);
    theta = zeros(1,n);
    psi = zeros(1,n);
    dtheta = zeros(1,n);
    dpsi = zeros(1,n);

    dr = zeros(3,1) ;
    rp = sc_states(1:3) ;
    r(:,1) = sc_states(1:3) ;
    v(:,1) = sc_states(4:6) ;
    t(:,1) = tspan(1) ;
    
%     % set orbit states with friendly names
%     a = states(1);
%     e = states(2);
%     i = states(3);
%     RAAN = states(4);
%     aop = states(5);
%     ta = states(6);

    % set attitude states with friendly names
    phi(1) = tether_state0(1);
    theta(1) = tether_state0(2);
    psi(1) = tether_state0(3);
    dphi(1) = tether_state0(4);
    dtheta(1) = tether_state0(5);
    dpsi(1) = tether_state0(6);
    
    % set tether parameters with friendly names
    L = tether_param(1);
    m1 = tether_param(2);
    m2 = tether_param(3);
    mt = tether_param(4);
    m = sum(tether_param(2:4));
    Ix = tether_param(5);
    Iy = tether_param(6);
    Iz = tether_param(7);
    current_type = tether_param(8);
    current_val = tether_param(9);
    
    ii = 2 ;
    dv = zeros(3,1) ;
    while t(ii-1) < tspan(2) && norm(r(:,ii-1)) >= (6378 + 100 ) 
%         
        [ r(:,ii) , v(:,ii) ] = NewStateUV( r(:,ii-1) , v(:,ii-1) , dt , mu  ) ; 
        eps = dot( r(:,ii) , dr )/norm( r(:,ii) )^2 ;
        if dr ~= zeros(3,1)
            f = ( 1/eps )*( 1 - ( 1/ ( 1 - 2*eps )^1.5 ) ) ;
        else
            f = 0 ;
        end

        if current_type == 0
            I = current_val;
        elseif current_type == 1
            I = OML(tether_param);
        elseif current_type == 2
            % while energy is less than required to reach limit from 
            limit_libration = 35;
            I = current_val;
            lim_e = m*L*1000*(1-cos(deg2rad(limit_libration)));
            pe_p = m*L*1000*(1-cos(phi(ii)));
            ke_p = (1/2)*Ix*(dphi(ii)*4)^2;
    %             ke_p = (1/2)*Ix*(dphi)^2;
            if phi < 0
                pe_p = -pe_p;
            end
            if dphi < 0
                ke_p = -ke_p;
            end
            e_p = pe_p + ke_p;
            pe_t = m*L*1000*(1-cos(theta));
            ke_t = (1/2)*Iy*(dtheta*4)^2;
    %             ke_t = (1/2)*Iy*(dtheta)^2;
            if theta < 0
                pe_t = -pe_t;
            end
            if dtheta < 0
                ke_t = -ke_t;
            end
            e_t = pe_t + ke_t;
            if abs(e_p) > lim_e
                I = 0;
            elseif abs(e_t) > lim_e
                I = 0;
            end
        end
    
        % Find the instantaneous Lorentz force, Lorentz torque, and gravity
        % gradient
        [Bx, By, Bz] = MagField_NonTilted(states);
    %         [Bx, By, Bz] = MagField_igrf(states, t);
        [fr, fs, fw] = edt_forces(states, tether_param, I, Bx, By, Bz);
        Tq = edt_torque(states, tether_param, I, Bx, By, Bz);
        Tgg = gravity_grad_torque(states, tether_param, mu);
        
        ap = [fr, fs, fw]'./m;
        da = ap*1e3 + ( mu/norm(r(:,ii))^3 )*( f*eps*rp - dr ) ;
        dv = da*dt + dv ;
        dr = .5*da*dt^2 + dv*dt + dr ;
            rp = r(:,ii) + dr ;
            vp = v(:,ii) + dv ;
            r(:,ii) = rp ;
            v(:,ii) = vp ;
            dr = zeros(3,1) ;
            dv = zeros(3,1) ;
        % attitude change from 16.1 of Spacecraft dynamics
        wo = -sqrt(mu/r^3);
        body_torque = Tgg+Tq;
        ddphi = dpsi*wo + ((Iz-Iy)*(wo^2*phi+wo*dpsi)+body_torque(1))/Ix;
        ddtheta = body_torque(2)/Iy;        
        ddpsi = 0;
        dstate = [da; de; di; dRAAN; daop; dta; dphi; dtheta; dpsi; ddphi; ddtheta; ddpsi];
        t(ii) = t(ii-1) + dt ;
        ii = ii + 1 ;
    end
    

  
    
    

    function [Bx, By, Bz] = MagField_NonTilted(states)
        % Earth’s magnetic field vector in the Euler–Hill frame whose 
        % components are defined by (assuming a nontilted dipole)

        % Define variables
        a = states(1);
        e = states(2);
        i = states(3);
        RAAN = states(4);
        aop = states(5);
        ta = states(6);
        p = a*(1-e^2);
        
        % Calculate magnetic field based on https://doi.org/10.2514/1.12016
        Bo = 3.12e-5;
        Bo_over_R3 = Bo/(6378/(p/(1+e*cos(ta))))^3;

        Bx = -2*Bo_over_R3*sin(aop + ta)*sin(i);
        By = Bo_over_R3*cos(aop+ta)*sin(i);
        Bz = Bo_over_R3*cos(i);
    end

    function [Bx, By, Bz] = MagField_igrf(states, t)
        % Earth’s magnetic field vector in the Euler–Hill frame whose 
        % components are defined by (assuming a nontilted dipole)

        % Define variables
        a = states(1);
        e = states(2);
        i = states(3);
        RAAN = states(4);
        aop = states(5);
        ta = states(6);
        p = a*(1-e^2);
        
        % calculate position a few ways
        [ rvec , vvec ] = coes2state( [sqrt(mu*a*(1-e^2)), i, e, RAAN, aop, ta] , mu );
        tsec = t;
        tmin = 0;
        thour = 12;
        tday = 1;
        while tsec >= 60
            tsec = tsec - 60;
            tmin = tmin + 1;
        end
        while tmin >= 60
            tmin = tmin - 60;
            thour = thour + 1;
        end
        while thour >= 24
            thour = thour - 24;
            tday = tday + 1;
        end
        time = [2020 1 tday thour tmin tsec];
        [ rvec_ecef , vvec_ecef ] = eci2ecef(time, rvec, vvec);
        LLA = ecef2lla(rvec_ecef'*1000);
%         LLA = eci2lla(rvec'.*1000,time);
        if LLA(3) < 0
            LLA(3) = (norm(rvec)-6378)*1000;
        end
        [B_NED,~,~,~,~] = igrfmagm(LLA(3),LLA(1),LLA(2),decyear(time),13);
        [B_ECEFX, B_ECEFY, B_ECEFZ] = ned2ecefv(B_NED(1),B_NED(2),B_NED(3),LLA(1),LLA(2));
        [r_eci, B_eci] = ecef2eci(time, rvec_ecef, [B_ECEFX, B_ECEFY, B_ECEFZ]);
        y_eh_ineci = r_eci/norm(r_eci);
        z_eh_ineci = cross(rvec,vvec)/norm(cross(rvec,vvec));
        x_eh_ineci = cross(y_eh_ineci,z_eh_ineci)/norm(cross(y_eh_ineci,z_eh_ineci));
        eci2eh = [x_eh_ineci, y_eh_ineci, z_eh_ineci];
        B_eh = eci2eh*B_eci;
        %NEED IN LVLH
        B_eh = B_eh*1e-9;
        Bx = B_eh(1);
        By = B_eh(2);
        Bz = B_eh(3);

        % Calculate magnetic field based on https://doi.org/10.2514/1.12016
%         Bo = 3.12e-5;
%         Bo_over_R3 = Bo/(6378/(p/(1+e*cos(ta))))^3;
% 
%         Bx = -2*Bo_over_R3*sin(aop + ta)*sin(i);
%         By = Bo_over_R3*cos(aop+ta)*sin(i);
%         Bz = Bo_over_R3*cos(i);
    end

    function [fr, ftheta, fh] = edt_forces(states, tether_param, I,  Bx, By, Bz)
        % assumes a uniform current flowing in the tether
        
        % set states with friendly names
        theta = states(8);
        phi = states(7);

        L = tether_param(1);
        m = sum(tether_param(2:4));

        

        fr = I*L/m*(Bz*sin(theta)*cos(phi)-By*sin(phi));
        ftheta = I*L/m*(Bx*sin(phi)-Bz*cos(theta)*cos(phi));
        fh = I*L/m*(By*cos(theta)*cos(phi)-Bx*sin(theta)*cos(phi));
    end

    function Tq = edt_torque(states, tether_param, I, Bx, By, Bz)
        % assumes a uniform current flowing in the tether. In the case 
        % where the mass distribution is perfectly symmetrical, then the 
        % net torque about the center of mass is zero. This is an ideal 
        % scenario that is unlikely to be achieved in practice. This effect
        % will still be present in bare-wire tethers, but the relationship
        % is not as straightforward because of the nonuniform distribution 
        % of the electric current.
        phi = states(7);
        theta = states(8);     

        % set states with friendly names
        L = tether_param(1);
        m1 = tether_param(2);
        m2 = tether_param(3);
        mt = tether_param(4);
        m = sum(tether_param(2:4));

        % nondimensional parameter that defines a measure of the average 
        % moment arm of the electromagnetic torque acting on the tether 
        % relative to the system center of mass
        PHI = (m1^2-m2^2+mt*(m1-m2))/m^2; 

        Qtheta = (I*(L*1000)^2/2)*PHI*cos(phi)*...
            (sin(phi)*(Bx*cos(theta)+By*sin(theta))-Bz*cos(phi));
        Qphi = -(I*(L*1000)^2/2)*PHI*(By*cos(theta)-Bx*sin(theta));
        Tq = [Qphi; Qtheta; 0];
    end

    function Tgg = gravity_grad_torque(states, tether_param, mu)
        % Finds the gravity gradient torque effecting the tether

        % set states with friendly names
        a = states(1);
        e = states(2);
        i = states(3);
        RAAN = states(4);
        aop = states(5);
        ta = states(6);
        phi = states(7);
        theta = states(8);
        psi = states(9);
    
        % Make inertia matrix
        Ix = tether_param(5);
        Iy = tether_param(6);
        Iz = tether_param(7);
        inertia = [Ix, 0, 0; 0, Iy, 0; 0, 0, Iz];
        
        % Calculate the position and velocity vector
        [ rvec , vvec ] = coes2state( [sqrt(mu*a*(1-e^2)), i, e, RAAN, aop, ta] , mu );
        r = norm(rvec);

        % Dynamics from "Spacecraft Dynamics and Control"
        Cbo = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
                sin(phi)*sin(theta)*cos(psi)-cos(theta)*sin(psi), ...
                sin(phi)*sin(theta)*sin(psi)+cos(phi)*sin(psi), ...
                sin(phi)*cos(theta);
                cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), ...
                cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), ...
                cos(phi)*cos(theta)];
        rb = Cbo*[0;0;-r];

        % transform inertia into eci
        rbcross = [0, -rb(3), rb(2);
                    rb(3), 0, -rb(1);
                    -rb(2), rb(2), 0];
        Tgg = 3*mu/r^5 * rbcross*inertia*rb;
    end
    
    function I = OML(tether_param)
        % Finds the OML current
        e = 1.60217663e-19 ; % charge of an electron
        me = 9.1093837e-31; % mass of an electron
        R = tether_param(10); %Tether radius (m)
        L = tether_param(1)*1000; % Tether Length (m)
        p = R*2*pi; % Tether cross sectional parameter (m)
        phiP = 100;%Cylindrical probe bias
        N0 = 0;% ambient electron density
        
        
        I = e*N0*(L*p/pi)*sqrt(2*e*phiP/me);


    end
    
end