clear ; close all ; clc ;
mu = 4904.8695 ;
Rcb = 1737 ;
n = 1e2 ; % vary time
m = 36  ; % vary arrival point
o = 36 ; % vary start point
d2r = (pi/180) ;
figure
% raanOriginal = [360,90,180,270;360,90,180,270;360,90,180,270;360,90,180,270;360,90,180,270;360,90,180,270;360,90,180,270;360,90,180,270;360,90,180,270;360,90,180,270] ;
% raanAllPlanes = [358.342402331547,89.3061682992730,179.143998809005,269.489912522383;356.576636310261,88.1391639669331,178.176135530858,268.259175629228;355.553832163779,86.3979182477023,177.633773360775,266.568930538507;354.430748650960,84.9684930802336,177.078135598480,265.442312262913;352.506376173485,83.9351962566774,176.074814147475,264.523524481739;351.036151222523,82.9059796450653,175.266915020007,263.316567467890;350.096139082184,81.1634923523413,174.793183863336,261.581040923586;348.662016378542,79.6758825933906,174.061806934162,260.324860516120;346.897722297406,79.0528669347316,173.143478655532,259.712622050327;345.390405826372,77.7443827919510,172.359390133487,258.178259591948] ;
% raanDelta = raanAllPlanes - raanOriginal;
% raanDeltaAvg = zeros(6,1) ;
% for ii = 1:length( raanOriginal ) 
%     raanDeltaAvg(ii) = mean(raanDelta(ii,:)) ;
% end

load( 'YoD_FC_7500.mat' ) ;
raanDriftAvg = zeros( 1 , 38 ) ;
for ii = 2:38
    raanDriftAvg(ii) = mean( raanP(ii,:) - raanP(1,:) ) ;
end

for plane = 1:4
    for d = 2:38
    % set arrival / nominal orbit
    incNom = 55 ;
%     incNom = incP(d,plane) ;
    aopNom = 0 ; 
    eccNom = 1e-9 ;
    hNom = hP(1,plane) ;
    thetaS = 0 ;
    thetaEnd = 360 ;
    thetaNom = linspace( thetaS, thetaEnd , m ) ;
    % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
    rNom = zeros( 3 , m ) ;
    vNom = zeros( 3 , m ) ;
    raanNom = raanP(1,plane) + raanDriftAvg(d) ;
%     raanNom = raanP(d,plane) ;
        for ii = 1:m
            [ rNom(:,ii) , vNom(:,ii) ] = coes2state( [ hNom , incNom(1)*d2r , eccNom , raanNom*d2r , aopNom(1)*d2r , thetaNom(ii)*d2r ] , mu ) ;
        end

        % set beginning / perturbed orbit
    %     incPer = [55.2818133233726;56.3257284620972;56.9635291085618;56.8906611133049;57.3164432916295;58.3723838156408;58.9685428999927;59.0041674165460;59.3826542638042;60.434174268500510] ;
        
        smaPer = (hP(d,plane)^2)/(mu*(1-eccP(d,plane)^2)) ;
        T = 2*pi*sqrt( smaPer^3 / mu ) ;
        % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
        [ rPer , vPer ] = coes2state( [ hP(d,plane) , incP(d,plane)*d2r , eccP(d,plane) , raanP(d,plane)*d2r , aopP(d,plane)*d2r , thetaP(d,plane)*d2r ] , mu ) ;

        % find min time
        tp = zeros( 1 , m ) ;
        tmin = zeros( m , 2 ) ;
        for ii = 1:m
            [ ~ , ~ , ~ , tmin(ii,:) , tp(ii) ] = Lambert_MinEnergy( rPer , rNom(:,ii) , mu ) ;
        end

        % set time of flight
        tf = zeros( n , m ) ;
        for ii = 1:m
            tf(:,ii) = linspace( tp(ii) , T*2 , n ) ;
        end

        % % time plot
        % figure
        % hold on
        % plot( thetaNom/d2r , tp )
        % plot( thetaNom/d2r , tmin(:,1) )
        % plot( thetaNom/d2r , tmin(:,2) )
        % hold off

        % use lamberts solver
        va = zeros( 3 , n , m ) ;
        vb = zeros( 3 , n , m ) ;
        deltav = zeros( n , m ) ;
        collision = zeros( n , m ) ;
        for ii = 1:n
            for jj = 1:m
                [ va(:,ii,jj) , vb(:,ii,jj) , deltav(ii,jj) , collision(ii,jj) ] = Targeting_FF( rPer , rNom(:,jj) , vPer , vNom(:,jj) , tf(ii,jj) , mu , Rcb ) ;
                if collision(ii,jj) == 1 
                    deltav(ii,jj) = 10 ;
                end
            end
        end

    %     % plot of delta v
    %     figure
    %     surf( tf*(1/(60*60)) , thetaNom , deltav , 'EdgeColor' , 'none' )
    %     xlabel( 'Time of Flight [hours]' )
    %     ylabel( 'Difference in True Anomaly [degrees]' )
    %     zlabel( 'Delta-V [km/s]' )
    %     axis([ min(tp)*(1/(60*60)) , 2*T*(1/(60*60)) , 0 , 360 , 0 , 10 ])


        [ min_perCol , rows ] = min( deltav ) ;
        [ mindv(d) , col ] = min( min_perCol ) ;
        row = rows( col ) ;

        tol = 1e-8 ;
        [ t , rv ] = TwoBody( [ 0 , 1.1*T ] , rNom(:,col) , vNom(:,col) , mu , tol ) ;
        [ tPer , rvPer ] = TwoBody( [ 0 , 1.1*T ] , rPer , vPer , mu , tol ) ;
        [ tTran , rvTran ] = TwoBody( [ 0 , tf( row , col ) ] , rPer , va(:,row,col) , mu , tol ) ;

%         figure
%         plot3( rv(:,1) , rv(:,2) , rv(:,3) , 'k' )
%         axis equal
%         hold on
%         plot3( rvPer(:,1) , rvPer(:,2) , rvPer(:,3) , 'b' )
%         plot3( rvTran(:,1) , rvTran(:,2) , rvTran(:,3) , 'r' )
%         plot3( rPer(1) , rPer(2) , rPer(3) , '*k' )
%         plot3( rNom(1,col) , rNom(2,col) , rNom(3,col) , '*b' )
%         legend( 'Perturbed Orbit' , 'Nominal Orbit' , 'Transfer Trajectory' )

        dir = [ "x" ; "y" ; "z" ] ;
        burn = va(:,row,col) - vPer ;
        burn2 = vNom(:,col) - vb(:,row,col) ;
        dispva = [ burn(1) , dir(1) , burn(2) , dir(2) , burn(3) , dir(3) ] ;
        dispvb = [ burn2(1) , dir(1) , burn2(2) , dir(2) , burn2(3) , dir(3) ] ;
        fprintf( 'The delta-v is %f m/s \n' , deltav(row,col)*1e3 )
        fprintf( 'The time of transfer is %f hours \n' , tf(row,col)*(1/(60*60)) )
%         fprintf( 'The initial burn is %f km/s in %s \n' , dispva )
%         fprintf( 'The final burn is %f km/s in %s \n' , dispvb )

    end

    for ii = 1:37 
        daysnum(ii) = ii*10 ;
        dvyear(ii) = mindv(d)*(365/daysnum(ii)) ;
    end

    
    plot( daysnum , dvyear )
    axis([0,370,0,1])
    hold on

end
xlabel( 'Time to Correction [days]' )
ylabel( 'Delta-V [km]' )
legend( 'Plane 1' , 'Plane 2' , 'Plane 3' , 'Plane 4' )

%% Functions
function [ r , v ] = coes2state( COES , mu )
% h , inc , ecc , RAAN , omega , theta , a , rp , ra 
    
h = COES(1) ;
inc = COES(2) ;
ecc = COES(3) ;
RAAN = COES(4) ;
omega = COES(5) ;
theta = COES(6) ;

    r_peri = (h^2/mu) * ( 1/( 1 + ecc*cos(theta) ) ) * [ cos( theta ) ; sin( theta ) ; 0 ] ;
    v_peri = (mu/h) * [ -sin( theta ) ; ecc+cos(theta) ; 0 ] ;

    Q(1,1) = -sin(RAAN)*cos(inc)*sin(omega) + cos(RAAN)*cos(omega) ;
    Q(1,2) = -sin(RAAN)*cos(inc)*cos(omega) - cos(RAAN)*sin(omega) ;
    Q(1,3) = sin(RAAN)*sin(inc) ;
    Q(2,1) = cos(RAAN)*cos(inc)*sin(omega) + sin(RAAN)*cos(omega) ;
    Q(2,2) = cos(RAAN)*cos(inc)*cos(omega) - sin(RAAN)*sin(omega) ;
    Q(2,3) = -cos(RAAN)*sin(inc) ;
    Q(3,1) = sin(inc)*sin(omega) ;
    Q(3,2) = sin(inc)*cos(omega) ;
    Q(3,3) = cos(inc) ;

    r = Q*r_peri ;
    v = Q*v_peri ;
end

function [ v1 , amin , emin , tmin , tp ] = Lambert_MinEnergy( r1 , r2 , mu )
% finds minimum energy, or fundamental ellipse, between two points
% r1,2 are position vectors, columns, in km
% v1 is column vectors in km/s
% amin is semimajor axis of min energy ellipse
% emin is eccentricity
% tmin is time of transit, short way then long way
% tp is parabolic time of transit

r1mag = norm( r1 ) ;
r2mag = norm( r2 ) ;
cos_ta = dot( r1 , r2 )/( r1mag*r2mag ) ;
c = sqrt( r1mag^2 + r2mag^2 - 2*r1mag*r2mag*cos_ta ) ;
s = ( r1mag + r2mag + c )/2 ;
amin = s/2 ;
pmin = ( r1mag*r2mag / c )*( 1 - cos_ta ) ;
emin = sqrt( 1 - ( 2*pmin )/s ) ;
beta = 2*asin( sqrt( ( s - c )/s ) ) ;
tmin(1) = sqrt( amin^3 / mu )*( pi - ( beta - sin( beta ) ) ) ;
tmin(2) = sqrt( amin^3 / mu )*( pi + ( beta - sin( beta ) ) ) ;
tp = ( 1/3 )*sqrt( 2/mu )*( s^1.5 - ( s - c )^1.5 ) ;
v1 = sqrt( mu*pmin )/( r1mag*r2mag*sin( acos( cos_ta ) ) )*( r2 - ( 1 - ( r2mag/pmin )*( 1 - cos_ta ) )*r1 ) ;
    
end

function [ va , vb , deltav , collision ] = Targeting_FF( rint , rtgt , vint , vtgt , dt , mu , Rcb )
% finds low deltav lamberts solution from a fixed point state to a fixed
% state. No propagation of target through transit time. Checks for
% collision with central body
    Tran_n = cross( rint , rtgt ) ;
    h_n = cross( rint , vint ) ;
    if dot( Tran_n , h_n ) > 0
        tm = 1 ;
    else 
        tm = -1 ;
    end
    [ va , vb ] = Lambert_Izzo( rint , rtgt , dt , tm , mu ) ;
    collision = HitCB( rint , rtgt , va , vb , mu , Rcb) ;
    deltav = abs( norm( va - vint ) + norm( vb - vtgt ) ) ;
    
    function [ collision ] = HitCB( rint , rtgt , va , vb , mu , Rcb )
        if dot( rint , va ) > 0 && dot( rtgt , vb )
            coes = state2coes( [ rint ; va ] , mu ) ;
            if coes(8) < Rcb
                collision = 1 ;
            else
                collision = 0 ;
            end
        else
            collision = 0 ;
        end
    end
end

function [ V1 , V2 ] = Lambert_Izzo( r1 , r2 , tf , tm , mu ) 
    [ V1r , V2r , ~ , ~ ] = Lambert_IzzoFun( r1' , r2' , tm*tf/(24*60*60) , 0 , mu ) ;
    V1 = V1r' ;
    V2 = V2r' ;
    function [V1,...
              V2, ...
              extremal_distances,...
              exitflag] = Lambert_IzzoFun(r1vec,...
                                        r2vec,...
                                        tf,...
                                        m,...
                                        muC) %#coder
    % original documentation:
    %{
     This routine implements a new algorithm that solves Lambert's problem. The
     algorithm has two major characteristics that makes it favorable to other
     existing ones.
     1) It describes the generic orbit solution of the boundary condition
     problem through the variable X=log(1+cos(alpha/2)). By doing so the
     graph of the time of flight become defined in the entire real axis and
     resembles a straight line. Convergence is granted within few iterations
     for all the possible geometries (except, of course, when the transfer
     angle is zero). When multiple revolutions are considered the variable is
     X=tan(cos(alpha/2)*pi/2).
     2) Once the orbit has been determined in the plane, this routine
     evaluates the velocity vectors at the two points in a way that is not
     singular for the transfer angle approaching to pi (Lagrange coefficient
     based methods are numerically not well suited for this purpose).
     As a result Lambert's problem is solved (with multiple revolutions
     being accounted for) with the same computational effort for all
     possible geometries. The case of near 180 transfers is also solved
     efficiently.
      We note here that even when the transfer angle is exactly equal to pi
     the algorithm does solve the problem in the plane (it finds X), but it
     is not able to evaluate the plane in which the orbit lies. A solution
     to this would be to provide the direction of the plane containing the
     transfer orbit from outside. This has not been implemented in this
     routine since such a direction would depend on which application the
     transfer is going to be used in.
     please report bugs to dario.izzo@esa.int
    %}
    % adjusted documentation:
    %{
     By default, the short-way solution is computed. The long way solution
     may be requested by giving a negative value to the corresponding
     time-of-flight [tf].
     For problems with |m| > 0, there are generally two solutions. By
     default, the right branch solution will be returned. The left branch
     may be requested by giving a negative value to the corresponding
     number of complete revolutions [m].
    %}
    % Authors
    % .-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.
    % Name       : Dr. Dario Izzo
    % E-mail     : dario.izzo@esa.int
    % Affiliation: ESA / Advanced Concepts Team (ACT)
    % Made more readible and optimized for speed by Rody P.S. Oldenhuis
    % Code available in MGA.M on   http://www.esa.int/gsp/ACT/inf/op/globopt.htm
    % last edited 12/Dec/2009
    % ADJUSTED FOR EML-COMPILATION 24/Dec/2009
        % initial values
        tol = 1e-14;    bad = false;     days = 86400;
        % work with non-dimensional units
        r1 = sqrt(r1vec*r1vec.');  r1vec = r1vec/r1;
        V = sqrt(muC/r1);          r2vec = r2vec/r1;
        T = r1/V;                  tf    = tf*days/T; % also transform to seconds
        % relevant geometry parameters (non dimensional)
        mr2vec = sqrt(r2vec*r2vec.');
        % make 100% sure it's in (-1 <= dth <= +1)
        dth = acos( max(-1, min(1, (r1vec*r2vec.')/mr2vec)) );
        % decide whether to use the left or right branch (for multi-revolution
        % problems), and the long- or short way
        leftbranch = sign(m);   longway = sign(tf);
        m = abs(m);             tf = abs(tf);
        if (longway < 0), dth = 2*pi - dth; end
        % derived quantities
        c      = sqrt(1 + mr2vec^2 - 2*mr2vec*cos(dth)); % non-dimensional chord
        s      = (1 + mr2vec + c)/2;                     % non-dimensional semi-perimeter
        a_min  = s/2;                                    % minimum energy ellipse semi major axis
        Lambda = sqrt(mr2vec)*cos(dth/2)/s;              % lambda parameter (from BATTIN's book)
        crossprd = [r1vec(2)*r2vec(3) - r1vec(3)*r2vec(2),...
                    r1vec(3)*r2vec(1) - r1vec(1)*r2vec(3),...% non-dimensional normal vectors
                    r1vec(1)*r2vec(2) - r1vec(2)*r2vec(1)];
        mcr       = sqrt(crossprd*crossprd.');           % magnitues thereof
        nrmunit   = crossprd/mcr;                        % unit vector thereof
        % Initial values
        % ---------------------------------------------------------
        % ELMEX requires this variable to be declared OUTSIDE the IF-statement
        logt = log(tf); % avoid re-computing the same value
        % single revolution (1 solution)
        if (m == 0)
            % initial values
            inn1 = -0.5233;      % first initial guess
            inn2 = +0.5233;      % second initial guess
            x1   = log(1 + inn1);% transformed first initial guess
            x2   = log(1 + inn2);% transformed first second guess
            % multiple revolutions (0, 1 or 2 solutions)
            % the returned soltuion depends on the sign of [m]
        else
            % select initial values
            if (leftbranch < 0)
                inn1 = -0.5234; % first initial guess, left branch
                inn2 = -0.2234; % second initial guess, left branch
            else
                inn1 = +0.7234; % first initial guess, right branch
                inn2 = +0.5234; % second initial guess, right branch
            end
            x1 = tan(inn1*pi/2);% transformed first initial guess
            x2 = tan(inn2*pi/2);% transformed first second guess
        end
        % since (inn1, inn2) < 0, initial estimate is always ellipse
        xx   = [inn1, inn2];  aa = a_min./(1 - xx.^2);
        bbeta = longway * 2*asin(sqrt((s-c)/2./aa));
        % make 100.4% sure it's in (-1 <= xx <= +1)
        aalfa = 2*acos(  max(-1, min(1, xx)) );
        % evaluate the time of flight via Lagrange expression
        y12  = aa.*sqrt(aa).*((aalfa - sin(aalfa)) - (bbeta-sin(bbeta)) + 2*pi*m);
        % initial estimates for y
        if m == 0
            y1 = log(y12(1)) - logt;
            y2 = log(y12(2)) - logt;
        else
            y1 = y12(1) - tf;
            y2 = y12(2) - tf;
        end
        % Solve for x
        % ---------------------------------------------------------
        % Newton-Raphson iterations
        % NOTE - the number of iterations will go to infinity in case
        % m > 0  and there is no solution. Start the other routine in
        % that case
        err = inf;  iterations = 0; xnew = 0;
        while (err > tol)
            % increment number of iterations
            iterations = iterations + 1;
            % new x
            xnew = (x1*y2 - y1*x2) / (y2-y1);
            % copy-pasted code (for performance)
            if m == 0, x = exp(xnew) - 1; else x = atan(xnew)*2/pi; end
            a = a_min/(1 - x^2);
            if (x < 1) % ellipse
                beta = longway * 2*asin(sqrt((s-c)/2/a));
                % make 100.4% sure it's in (-1 <= xx <= +1)
                alfa = 2*acos( max(-1, min(1, x)) );
            else % hyperbola
                alfa = 2*acosh(x);
                beta = longway * 2*asinh(sqrt((s-c)/(-2*a)));
            end
            % evaluate the time of flight via Lagrange expression
            if (a > 0)
                tof = a*sqrt(a)*((alfa - sin(alfa)) - (beta-sin(beta)) + 2*pi*m);
            else
                tof = -a*sqrt(-a)*((sinh(alfa) - alfa) - (sinh(beta) - beta));
            end
            % new value of y
            if m ==0, ynew = log(tof) - logt; else ynew = tof - tf; end
            % save previous and current values for the next iterarion
            % (prevents getting stuck between two values)
            x1 = x2;  x2 = xnew;
            y1 = y2;  y2 = ynew;
            % update error
            err = abs(x1 - xnew);
            % escape clause
            if (iterations > 15), bad = true; break; end
        end
        % If the Newton-Raphson scheme failed, try to solve the problem
        % with the other Lambert targeter.
        if bad
            % NOTE: use the original, UN-normalized quantities
            [V1, V2, extremal_distances, exitflag] = ...
                lambert_LancasterBlanchard(r1vec*r1, r2vec*r1, longway*tf*T, leftbranch*m, muC);
            return
        end
        % convert converged value of x
        if m==0, x = exp(xnew) - 1; else x = atan(xnew)*2/pi; end
        %{
          The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
          now need the conic. As for transfer angles near to pi the Lagrange-
          coefficients technique goes singular (dg approaches a zero/zero that is
          numerically bad) we here use a different technique for those cases. When
          the transfer angle is exactly equal to pi, then the ih unit vector is not
          determined. The remaining equations, though, are still valid.
        %}
        % Solution for the semi-major axis
        a = a_min/(1-x^2);
        % Calculate psi
        if (x < 1) % ellipse
            beta = longway * 2*asin(sqrt((s-c)/2/a));
            % make 100.4% sure it's in (-1 <= xx <= +1)
            alfa = 2*acos( max(-1, min(1, x)) );
            psi  = (alfa-beta)/2;
            eta2 = 2*a*sin(psi)^2/s;
            eta  = sqrt(eta2);
        else       % hyperbola
            beta = longway * 2*asinh(sqrt((c-s)/2/a));
            alfa = 2*acosh(x);
            psi  = (alfa-beta)/2;
            eta2 = -2*a*sinh(psi)^2/s;
            eta  = sqrt(eta2);
        end
        % unit of the normalized normal vector
        ih = longway * nrmunit;
        % unit vector for normalized [r2vec]
        r2n = r2vec/mr2vec;
        % cross-products
        % don't use cross() (emlmex() would try to compile it, and this way it
        % also does not create any additional overhead)
        crsprd1 = [ih(2)*r1vec(3)-ih(3)*r1vec(2),...
                   ih(3)*r1vec(1)-ih(1)*r1vec(3),...
                   ih(1)*r1vec(2)-ih(2)*r1vec(1)];
        crsprd2 = [ih(2)*r2n(3)-ih(3)*r2n(2),...
                   ih(3)*r2n(1)-ih(1)*r2n(3),...
                   ih(1)*r2n(2)-ih(2)*r2n(1)];
        % radial and tangential directions for departure velocity
        Vr1 = 1/eta/sqrt(a_min) * (2*Lambda*a_min - Lambda - x*eta);
        Vt1 = sqrt(mr2vec/a_min/eta2 * sin(dth/2)^2);
        % radial and tangential directions for arrival velocity
        Vt2 = Vt1/mr2vec;
        Vr2 = (Vt1 - Vt2)/tan(dth/2) - Vr1;
        % terminal velocities
        V1 = (Vr1*r1vec + Vt1*crsprd1)*V;
        V2 = (Vr2*r2n + Vt2*crsprd2)*V;
        % exitflag
        exitflag = 1; % (success)
        % also compute minimum distance to central body
        % NOTE: use un-transformed vectors again!
        extremal_distances = ...
            minmax_distances(r1vec*r1, r1, r2vec*r1, mr2vec*r1, dth, a*r1, V1, V2, m, muC);
    end
    % -----------------------------------------------------------------
    % Lancaster & Blanchard version, with improvements by Gooding
    % Very reliable, moderately fast for both simple and complicated cases
    % -----------------------------------------------------------------
    function [V1,...
              V2,...
              extremal_distances,...
              exitflag] = lambert_LancasterBlanchard(r1vec,...
                                                     r2vec,...
                                                     tf,...
                                                     m,...
                                                     muC) %#coder
    %{
    LAMBERT_LANCASTERBLANCHARD       High-Thrust Lambert-targeter
    lambert_LancasterBlanchard() uses the method developed by
    Lancaster & Blancard, as described in their 1969 paper. Initial values,
    and several details of the procedure, are provided by R.H. Gooding,
    as described in his 1990 paper.
    %}
    % Please report bugs and inquiries to:
    %
    % Name       : Rody P.S. Oldenhuis
    % E-mail     : oldenhuis@gmail.com
    % Licence    : 2-clause BSD (see License.txt)
    % If you find this work useful, please consider a donation:
    % https://www.paypal.me/RodyO/3.5
        % ADJUSTED FOR EML-COMPILATION 29/Sep/2009
        % manipulate input
        tol     = 1e-12;                            % optimum for numerical noise v.s. actual precision
        r1      = sqrt(r1vec*r1vec.');              % magnitude of r1vec
        r2      = sqrt(r2vec*r2vec.');              % magnitude of r2vec
        r1unit  = r1vec/r1;                         % unit vector of r1vec
        r2unit  = r2vec/r2;                         % unit vector of r2vec
        crsprod = cross(r1vec, r2vec, 2);           % cross product of r1vec and r2vec
        mcrsprd = sqrt(crsprod*crsprod.');          % magnitude of that cross product
        th1unit = cross(crsprod/mcrsprd, r1unit);   % unit vectors in the tangential-directions
        th2unit = cross(crsprod/mcrsprd, r2unit);
        % make 100.4% sure it's in (-1 <= x <= +1)
        dth = acos( max(-1, min(1, (r1vec*r2vec.')/r1/r2)) ); % turn angle
        % if the long way was selected, the turn-angle must be negative
        % to take care of the direction of final velocity
        longway = sign(tf); tf = abs(tf);
        if (longway < 0), dth = dth-2*pi; end
        % left-branch
        leftbranch = sign(m); m = abs(m);
        % define constants
        c  = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(dth));
        s  = (r1 + r2 + c) / 2;
        T  = sqrt(8*muC/s^3) * tf;
        q  = sqrt(r1*r2)/s * cos(dth/2);
        % general formulae for the initial values (Gooding)
        % -------------------------------------------------
        % some initial values
        T0  = LancasterBlanchard(0, q, m);
        Td  = T0 - T;
        phr = mod(2*atan2(1 - q^2, 2*q), 2*pi);
        % initial output is pessimistic
        V1 = NaN(1,3);    V2 = V1;    extremal_distances = [NaN, NaN];
        % single-revolution case
        if (m == 0)
            x01 = T0*Td/4/T;
            if (Td > 0)
                x0 = x01;
            else
                x01 = Td/(4 - Td);
                x02 = -sqrt( -Td/(T+T0/2) );
                W   = x01 + 1.7*sqrt(2 - phr/pi);
                if (W >= 0)
                    x03 = x01;
                else
                    x03 = x01 + (-W).^(1/16).*(x02 - x01);
                end
                lambda = 1 + x03*(1 + x01)/2 - 0.03*x03^2*sqrt(1 + x01);
                x0 = lambda*x03;
            end
            % this estimate might not give a solution
            if (x0 < -1), exitflag = -1; return; end
        % multi-revolution case
        else
            % determine minimum Tp(x)
            xMpi = 4/(3*pi*(2*m + 1));
            if (phr < pi)
                xM0 = xMpi*(phr/pi)^(1/8);
            elseif (phr > pi)
                xM0 = xMpi*(2 - (2 - phr/pi)^(1/8));
            % EMLMEX requires this one
            else
                xM0 = 0;
            end
            % use Halley's method
            xM = xM0;  Tp = inf;  iterations = 0;
            while abs(Tp) > tol
                % iterations
                iterations = iterations + 1;
                % compute first three derivatives
                [dummy, Tp, Tpp, Tppp] = LancasterBlanchard(xM, q, m);%#ok
                % new value of xM
                xMp = xM;
                xM  = xM - 2*Tp.*Tpp ./ (2*Tpp.^2 - Tp.*Tppp);
                % escape clause
                if mod(iterations, 7), xM = (xMp+xM)/2; end
                % the method might fail. Exit in that case
                if (iterations > 25), exitflag = -2; return; end
            end
            % xM should be elliptic (-1 < x < 1)
            % (this should be impossible to go wrong)
            if (xM < -1) || (xM > 1), exitflag = -1; return; end
            % corresponding time
            TM = LancasterBlanchard(xM, q, m);
            % T should lie above the minimum T
            if (TM > T), exitflag = -1; return; end
            % find two initial values for second solution (again with lambda-type patch)
            % --------------------------------------------------------------------------
            % some initial values
            TmTM  = T - TM;   T0mTM = T0 - TM;
            [dummy, Tp, Tpp] = LancasterBlanchard(xM, q, m);%#ok
            % first estimate (only if m > 0)
            if leftbranch > 0
                x   = sqrt( TmTM / (Tpp/2 + TmTM/(1-xM)^2) );
                W   = xM + x;
                W   = 4*W/(4 + TmTM) + (1 - W)^2;
                x0  = x*(1 - (1 + m + (dth - 1/2)) / ...
                    (1 + 0.15*m)*x*(W/2 + 0.03*x*sqrt(W))) + xM;
                % first estimate might not be able to yield possible solution
                if (x0 > 1), exitflag = -1; return; end
            % second estimate (only if m > 0)
            else
                if (Td > 0)
                    x0 = xM - sqrt(TM/(Tpp/2 - TmTM*(Tpp/2/T0mTM - 1/xM^2)));
                else
                    x00 = Td / (4 - Td);
                    W = x00 + 1.7*sqrt(2*(1 - phr));
                    if (W >= 0)
                        x03 = x00;
                    else
                        x03 = x00 - sqrt((-W)^(1/8))*(x00 + sqrt(-Td/(1.5*T0 - Td)));
                    end
                    W      = 4/(4 - Td);
                    lambda = (1 + (1 + m + 0.24*(dth - 1/2)) / ...
                        (1 + 0.15*m)*x03*(W/2 - 0.03*x03*sqrt(W)));
                    x0     = x03*lambda;
                end
                % estimate might not give solutions
                if (x0 < -1), exitflag = -1; return; end
            end
        end
        % find root of Lancaster & Blancard's function
        % --------------------------------------------
        % (Halley's method)
        x = x0; Tx = inf; iterations = 0;
        while abs(Tx) > tol
            % iterations
            iterations = iterations + 1;
            % compute function value, and first two derivatives
            [Tx, Tp, Tpp] = LancasterBlanchard(x, q, m);
            % find the root of the *difference* between the
            % function value [T_x] and the required time [T]
            Tx = Tx - T;
            % new value of x
            xp = x;
            x  = x - 2*Tx*Tp ./ (2*Tp^2 - Tx*Tpp);
            % escape clause
            if mod(iterations, 7), x = (xp+x)/2; end
            % Halley's method might fail
            if iterations > 25, exitflag = -2; return; end
        end
        % calculate terminal velocities
        % -----------------------------
        % constants required for this calculation
        gamma = sqrt(muC*s/2);
        if (c == 0)
            sigma = 1;
            rho   = 0;
            z     = abs(x);
        else
            sigma = 2*sqrt(r1*r2/(c^2)) * sin(dth/2);
            rho   = (r1 - r2)/c;
            z     = sqrt(1 + q^2*(x^2 - 1));
        end
        % radial component
        Vr1    = +gamma*((q*z - x) - rho*(q*z + x)) / r1;
        Vr1vec = Vr1*r1unit;
        Vr2    = -gamma*((q*z - x) + rho*(q*z + x)) / r2;
        Vr2vec = Vr2*r2unit;
        % tangential component
        Vtan1    = sigma * gamma * (z + q*x) / r1;
        Vtan1vec = Vtan1 * th1unit;
        Vtan2    = sigma * gamma * (z + q*x) / r2;
        Vtan2vec = Vtan2 * th2unit;
        % Cartesian velocity
        V1 = Vtan1vec + Vr1vec;
        V2 = Vtan2vec + Vr2vec;
        % exitflag
        exitflag = 1; % (success)
        % also determine minimum/maximum distance
        a = s/2/(1 - x^2); % semi-major axis
        extremal_distances = minmax_distances(r1vec, r1, r1vec, r2, dth, a, V1, V2, m, muC);
    end
    % Lancaster & Blanchard's function, and three derivatives thereof
    function [T, Tp, Tpp, Tppp] = LancasterBlanchard(x, q, m)
        % protection against idiotic input
        if (x < -1) % impossible; negative eccentricity
            x = abs(x) - 2;
        elseif (x == -1) % impossible; offset x slightly
            x = x + eps;
        end
        % compute parameter E
        E  = x*x - 1;
        % T(x), T'(x), T''(x)
        if x == 1 % exactly parabolic; solutions known exactly
            % T(x)
            T = 4/3*(1-q^3);
            % T'(x)
            Tp = 4/5*(q^5 - 1);
            % T''(x)
            Tpp = Tp + 120/70*(1 - q^7);
            % T'''(x)
            Tppp = 3*(Tpp - Tp) + 2400/1080*(q^9 - 1);
        elseif abs(x-1) < 1e-2 % near-parabolic; compute with series
            % evaluate sigma
            [sig1, dsigdx1, d2sigdx21, d3sigdx31] = sigmax(-E);
            [sig2, dsigdx2, d2sigdx22, d3sigdx32] = sigmax(-E*q*q);
            % T(x)
            T = sig1 - q^3*sig2;
            % T'(x)
            Tp = 2*x*(q^5*dsigdx2 - dsigdx1);
            % T''(x)
            Tpp = Tp/x + 4*x^2*(d2sigdx21 - q^7*d2sigdx22);
            % T'''(x)
            Tppp = 3*(Tpp-Tp/x)/x + 8*x*x*(q^9*d3sigdx32 - d3sigdx31);
        else % all other cases
            % compute all substitution functions
            y  = sqrt(abs(E));
            z  = sqrt(1 + q^2*E);
            f  = y*(z - q*x);
            g  = x*z - q*E;
            % BUGFIX: (Simon Tardivel) this line is incorrect for E==0 and f+g==0
            % d  = (E < 0)*(atan2(f, g) + pi*m) + (E > 0)*log( max(0, f + g) );
            % it should be written out like so:
            if (E<0)
                d = atan2(f, g) + pi*m;
            elseif (E==0)
                d = 0;
            else
                d = log(max(0, f+g));
            end
            % T(x)
            T = 2*(x - q*z - d/y)/E;
            %  T'(x)
            Tp = (4 - 4*q^3*x/z - 3*x*T)/E;
            % T''(x)
            Tpp = (-4*q^3/z * (1 - q^2*x^2/z^2) - 3*T - 3*x*Tp)/E;
            % T'''(x)
            Tppp = (4*q^3/z^2*((1 - q^2*x^2/z^2) + 2*q^2*x/z^2*(z - x)) - 8*Tp - 7*x*Tpp)/E;
        end
    end
    % series approximation to T(x) and its derivatives
    % (used for near-parabolic cases)
    function [sig, dsigdx, d2sigdx2, d3sigdx3] = sigmax(y)
        % preload the factors [an]
        % (25 factors is more than enough for 16-digit accuracy)
        persistent an;
        if isempty(an)
            an = [
                4.000000000000000e-001;     2.142857142857143e-001;     4.629629629629630e-002
                6.628787878787879e-003;     7.211538461538461e-004;     6.365740740740740e-005
                4.741479925303455e-006;     3.059406328320802e-007;     1.742836409255060e-008
                8.892477331109578e-010;     4.110111531986532e-011;     1.736709384841458e-012
                6.759767240041426e-014;     2.439123386614026e-015;     8.203411614538007e-017
                2.583771576869575e-018;     7.652331327976716e-020;     2.138860629743989e-021
                5.659959451165552e-023;     1.422104833817366e-024;     3.401398483272306e-026
                7.762544304774155e-028;     1.693916882090479e-029;     3.541295006766860e-031
                7.105336187804402e-033];
        end
        % powers of y
        powers = y.^(1:25);
        % sigma itself
        sig = 4/3 + powers*an;
        % dsigma / dx (derivative)
        dsigdx = ( (1:25).*[1, powers(1:24)] ) * an;
        % d2sigma / dx2 (second derivative)
        d2sigdx2 = ( (1:25).*(0:24).*[1/y, 1, powers(1:23)] ) * an;
        % d3sigma / dx3 (third derivative)
        d3sigdx3 = ( (1:25).*(0:24).*(-1:23).*[1/y/y, 1/y, 1, powers(1:22)] ) * an;
    end
    % -----------------------------------------------------------------
    % Helper functions
    % -----------------------------------------------------------------
    % compute minimum and maximum distances to the central body
    function extremal_distances = minmax_distances(r1vec, r1,...
                                                   r2vec, r2,...
                                                   dth,...
                                                   a,...
                                                   V1, V2,...
                                                   m,...
                                                   muC)
        % default - minimum/maximum of r1,r2
        minimum_distance = min(r1,r2);
        maximum_distance = max(r1,r2);
        % was the longway used or not?
        longway = abs(dth) > pi;
        % eccentricity vector (use triple product identity)
        evec = ((V1*V1.')*r1vec - (V1*r1vec.')*V1)/muC - r1vec/r1;
        % eccentricity
        e = sqrt(evec*evec.');
        % apses
        pericenter = a*(1-e);
        apocenter  = inf;                    % parabolic/hyperbolic case
        if (e < 1), apocenter = a*(1+e); end % elliptic case
        % since we have the eccentricity vector, we know exactly where the
        % pericenter lies. Use this fact, and the given value of [dth], to
        % cross-check if the trajectory goes past it
        if (m > 0) % obvious case (always elliptical and both apses are traversed)
            minimum_distance = pericenter;
            maximum_distance = apocenter;
        else % less obvious case
            % compute theta1&2 ( use (AxB)-(CxD) = (C·B)(D·A) - (C·A)(B·D) ))
            pm1 = sign( r1*r1*(evec*V1.') - (r1vec*evec.')*(r1vec*V1.') );
            pm2 = sign( r2*r2*(evec*V2.') - (r2vec*evec.')*(r2vec*V2.') );
            % make 100.4% sure it's in (-1 <= theta12 <= +1)
            theta1 = pm1*acos( max(-1, min(1, (r1vec/r1)*(evec/e).')) );
            theta2 = pm2*acos( max(-1, min(1, (r2vec/r2)*(evec/e).')) );
            % points 1&2 are on opposite sides of the symmetry axis -- minimum
            % and maximum distance depends both on the value of [dth], and both
            % [theta1] and [theta2]
            if (theta1*theta2 < 0)
                % if |th1| + |th2| = turnangle, we know that the pericenter was
                % passed
                if abs(abs(theta1) + abs(theta2) - dth) < 5*eps(dth)
                    minimum_distance = pericenter;
                % this condition can only be false for elliptic cases, and
                % when it is indeed false, we know that the orbit passed
                % apocenter
                else
                    maximum_distance = apocenter;
                end
            % points 1&2 are on the same side of the symmetry axis. Only if the
            % long-way was used are the min. and max. distances different from
            % the min. and max. values of the radii (namely, equal to the apses)
            elseif longway
                minimum_distance = pericenter;
                if (e < 1), maximum_distance = apocenter; end
            end
        end
        % output argument
        extremal_distances = [minimum_distance, maximum_distance];
    end
end

function [ t , r , v ] = TwoBody( tspan , r0 , v0 , mu , tol )
    opts = odeset( 'RelTol' , tol , 'AbsTol' , tol ) ;
    [ t , rv ] = ode45( @TwoBodyForceFun , tspan , [ r0 ; v0 ] , opts , mu ) ;
    r = rv(:,1:3)' ;
    v = rv(:,4:6)' ;
    
    function drdv = TwoBodyForceFun( t , rv , mu )
        rvec = rv(1:3) ;
        vvec = rv(4:6) ;

        dr = vvec ;
        dv = ( -mu / norm(rvec)^3 )*rvec ;
        drdv = [ dr ; dv ] ;
    end

end

function COES = state2coes( state, mu )
% all angles output in radians
% h , inc , ecc , RAAN , omega , theta , a , rp , ra 
    
    Kh = [ 0 0 1 ] ; % K hat
    R = state(1:3) ;
    V = state(4:6) ;
    r = norm( R ) ;
    v = norm( V ) ;
    vr = dot( R , V )/r ; % radial velocity
    H = cross( R , V ) ; % specific angular momentum
    h = norm( H ) ; % specific angular momentum
    inc = acos( H(3)/h ) ; %inclination
    ECC = (1/mu)*( ( v^2 - mu/r )*R - r*vr*V ) ; %eccentricity vector
    ecc = norm( ECC ) ; % eccentricity
    N = cross( Kh , H ) ; % Node line
    n = norm( N ) ;
if n ~= 0
    RAAN = acos(N(1)/n) ; %Right ascension of ascending node
    if N(2) < 0 
        RAAN = 2*pi - RAAN ; %Right ascension of ascending node
    end
else
    RAAN = 0 ;
end 
    
if n ~= 0 
    if ecc >= 0 
        omega = acos(dot(N,ECC)/(n*ecc)) ; % Argument of perigee
        if ECC(3) < 0 
            omega = 2*pi - omega ; % Argument of perigee
        end
    else
        omega = 0 ;
    end
else
    omega = 0 ;
end
    
if ecc > 0
    theta = acos( dot( ECC , R )/( ecc*r ) ) ;         
    if vr < 0 
        theta = 2*pi - theta ; 
    end
else
    cp = cross( N , R ) ;
    if cp(3) >= 0
        theta = acos( dot( N , R )/( n*r ) ) ;
    else
        theta = 2*pi - acos( dot( N , R )/( n*r ) ) ;
    end
end

    a = (h^2)/( mu*( 1 - ecc^2 ) ) ; % semi-major axis
    rp = a*( 1 - ecc ) ;
    ra = a*( 1 + ecc ) ;

    
    COES = [ h , inc , ecc , RAAN , omega , theta , a , rp , ra ] ;
end

