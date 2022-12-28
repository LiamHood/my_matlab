function [ t , states] = BasicTether( tspan , sc_state0, tether_state0, tether_param, mu , tol )
% Uses Classical Orbital Elements a (semimajor axis), e (eccentricity), 
% i (inclination), RAAN (right ascension of ascending node), aop (argument
% of periapsis), ta (true anomaly)
%
% Considers electrodynamic force and gravity gradient
% will at some point maybe: atmospheric drag, spherical harmonics,
% third body of sun and moon, srp
% tether is treated as dumbbell model, rigid with lumped masses at the end

% tether_param = L, m1, m2, mt

    opts = odeset('RelTol', tol, 'AbsTol', tol ) ;
    [ t , states ] = ode45(@gauss_variations, tspan , ...
        [ sc_state0 ; tether_state0], opts, tether_param, mu) ;


    function dstate = gauss_variations(t, states, tether_param, mu)
        a = states(1);
        e = states(2);
        i = states(3);
        RAAN = states(4);
        aop = states(5);
        ta = states(6);

        phi = states(7);
        theta = states(8);
        psi = states(9);
        dphi = states(10);
        dtheta = states(11);
        dpsi = states(12);
        
        L = tether_param(1);
        m1 = tether_param(2);
        m2 = tether_param(3);
        mt = tether_param(4);
        m = sum(tether_param(2:4));
        Ix = tether_param(5);
        Iy = tether_param(6);
        Iz = tether_param(7);

        [fr, fs, fw] = edt_forces(states, tether_param);
        Tq = edt_torque(states, tether_param);
        Tgg = gravity_grad_torque(states, tether_param, mu);

        p = a*(1-e^2);
        h = sqrt(mu*p);
        r = norm(states(1:3));
        n = sqrt(mu/a^3);
        u = aop + ta;

        % orbital elements using Vallado p.636
        da = (2/(n*sqrt(1-e^2)))*(e*sin(ta)*fr + (p/r)*fs);
        de = (sqrt(1-e^2)/(n*a))*(sin(ta)*fr + ...
            (cos(ta)+(e+cos(ta))/(1+e*cos(ta)))*fs);
        di = (r*cos(u)/(n*a^2*sqrt(1-e^2)))*fw;
        dRAAN = (r*sin(u)/(n*a^2*sqrt(1-e^2)*sin(i)))*fw;
        daop = (sqrt(1-e^2)/(n*a*e))*(-cos(ta)*fr + ...
            sin(ta)*(1+r/p)*fs) - (r*cot(i)*sin(u)/h)*fw;
        dta = h/r^2 + (1/(e*h))*(p*cos(ta)*fr - (p+r)*sin(ta)*fs);
        
        % dynamics from 16.1 of Spacecraft dynamics
        wo = -sqrt(mu/r^3);
        body_torque = Tgg+Tq;
        ddphi = dpsi*wo + ((Iz-Iy)*(wo^2*phi+wo*dpsi)+body_torque(1))/Ix;
        ddtheta = body_torque(2)/Iy;        
        ddpsi = 0;
        dstate = [da; de; di; dRAAN; daop; dta; dphi; dtheta; dpsi; ddphi; ddtheta; ddpsi];
    end

    function [Bx, By, Bz] = MagField_NonTilted(states)
        % Earth’s magnetic field vector in the Euler–Hill frame whose 
        % components are defined by (assuming a nontilted dipole)
        a = states(1);
        e = states(2);
        i = states(3);
        RAAN = states(4);
        aop = states(5);
        ta = states(6);
        p = a*(1-e^2);
        
        Bo = 3.12e-5;
        Bo_over_R3 = Bo/(6378/(p/(1+e*cos(ta))))^3;

        Bx = -2*Bo_over_R3*sin(aop + ta)*sin(i);
        By = Bo_over_R3*cos(aop+ta)*sin(i);
        Bz = Bo_over_R3*cos(i);
    end

    function [fr, ftheta, fh] = edt_forces(states, tether_param)
        theta = states(8);
        phi = states(7);
        current_type = tether_param(8);
        if current_type == 0
            I = tether_param(9);
        elseif current_type == 1
            I = tether_param(9);
        end

        L = tether_param(1);
        m = sum(tether_param(2:4));

        [Bx, By, Bz] = MagField_NonTilted(states);

        fr = I*L/m*(Bz*sin(theta)*cos(phi)-By*sin(phi));
        ftheta = I*L/m*(Bx*sin(phi)-Bz*cos(theta)*cos(phi));
        fh = I*L/m*(By*cos(theta)*cos(phi)-Bx*sin(theta)*cos(phi));
    end

    function Tq = edt_torque(states, tether_param)
        % assumes a uniform current flowing in the tether. In the case 
        % where the mass distribution is perfectly symmetrical, then the 
        % net torque about the center of mass is zero. This is an ideal 
        % scenario that is unlikely to be achieved in practice. This effect
        % will still be present in bare-wire tethers, but the relationship
        % is not as straightforward because of the nonuniform distribution 
        % of the electric current.
        phi = states(7);
        theta = states(8);     
        current_type = tether_param(8);
        if current_type == 0
            I = tether_param(9);
        elseif current_type == 1
            I = tether_param(9);
        end

        L = tether_param(1);
        m1 = tether_param(2);
        m2 = tether_param(3);
        mt = tether_param(4);
        m = sum(tether_param(2:4));

        % nondimensional parameter that defines a measure of the average 
        % moment arm of the electromagnetic torque acting on the tether 
        % relative to the system center of mass
        PHI = (m1^2-m2^2+mt*(m1-m2))/m^2; 
        [Bx, By, Bz] = MagField_NonTilted(states);

        Qtheta = (I*(L*1000)^2/2)*PHI*cos(phi)*...
            (sin(phi)*(Bx*cos(theta)+By*sin(theta))-Bz*cos(phi));
        Qphi = -(I*(L*1000)^2/2)*PHI*(By*cos(theta)-Bx*sin(theta));
        Tq = [Qphi; Qtheta; 0];
    end

    function Tgg = gravity_grad_torque(states, tether_param, mu)
        a = states(1);
        e = states(2);
        i = states(3);
        RAAN = states(4);
        aop = states(5);
        ta = states(6);
        phi = states(7);
        theta = states(8);
        psi = states(9);

        Ix = tether_param(5);
        Iy = tether_param(6);
        Iz = tether_param(7);
        inertia = [Ix, 0, 0; 0, Iy, 0; 0, 0, Iz];
        
        [ rvec , vvec ] = coes2state( [sqrt(mu*a*(1-e^2)), i, e, RAAN, aop, ta] , mu );
        r = norm(rvec);

%         Dynamics from "Spacecraft Dynamics and Control"
        Cbo = [cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta);
                sin(phi)*sin(theta)*cos(psi)-cos(theta)*sin(psi), ...
                sin(phi)*sin(theta)*sin(psi)+cos(phi)*sin(psi), ...
                sin(phi)*cos(theta);
                cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi), ...
                cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi), ...
                cos(phi)*cos(theta)];
        rb = Cbo*[0;0;-r];

        % need L in terms ECI to get r new
        rbcross = [0, -rb(3), rb(2);
                    rb(3), 0, -rb(1);
                    -rb(2), rb(2), 0];
        Tgg = 3*mu/r^5 * rbcross*inertia*rb;
    end
    
    function I = OML(tether_param)
        R = tether_param(10); %Tether radius (m)
        L = tether_param(1)*1000; % Tether Length (m)
        phiP = 100;%Cylindrical probe bias
        N0 = 0;% ambient electron density
        
        
        Ith = 2*pi*R*L*e*N0*sqrt(k)*Te/(2*pi*me);
        I = Ith*sqrt(4*e*phiP/(pi*k*Te));

    end
    
end