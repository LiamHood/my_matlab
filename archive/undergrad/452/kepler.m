    % ------------------------------------------------------------------------------
    %
    %                           function kepler
    %
    %  this function solves keplers problem for orbit determination and returns a
    %    future geocentric equatorial (ijk) position and velocity vector.  the
    %    solution uses universal variables.
    %
    %  author        : david vallado                  719-573-2600   22 jun 2002
    %
    %  revisions
    %    vallado     - fix some mistakes                             13 apr 2004
    %
    %  inputs          description                    range / units
    %    ro          - ijk position vector - initial  km
    %    vo          - ijk velocity vector - initial  km / s
    %    dtsec       - length of time to propagate    s
    %
    %  outputs       :
    %    r           - ijk position vector            km
    %    v           - ijk velocity vector            km / s
    %    error       - error flag                     'ok', ...
    %
    %  locals        :
    %    f           - f expression
    %    g           - g expression
    %    fdot        - f dot expression
    %    gdot        - g dot expression
    %    xold        - old universal variable x
    %    xoldsqrd    - xold squared
    %    xnew        - new universal variable x
    %    xnewsqrd    - xnew squared
    %    znew        - new value of z
    %    c2new       - c2(psi) function
    %    c3new       - c3(psi) function
    %    dtsec       - change in time                 s
    %    timenew     - new time                       s
    %    rdotv       - result of ro dot vo
    %    a           - semi or axis                   km
    %    alpha       - reciprocol  1/a
    %    sme         - specific mech energy           km2 / s2
    %    period      - time period for satellite      s
    %    s           - variable for parabolic case
    %    w           - variable for parabolic case
    %    h           - angular momentum vector
    %    temp        - temporary real*8 value
    %    i           - index
    %
    %  coupling      :
    %    mag         - magnitude of a vector
    %    findc2c3    - find c2 and c3 functions
    %
    %  references    :
    %    vallado       2004, 95-103, alg 8, ex 2-4
    %
    % [r, v] =  kepler  ( ro, vo, dtsec );
    % ------------------------------------------------------------------------------

    function [r, v] =  kepler  ( ro, vo, dtseco )
    %function [r,v,errork] =  kepler  ( ro,vo, dtseco, fid );

    % -------------------------  implementation   -----------------
    % set constants and intermediate printouts
mu = 398600.4418 ;
    numiter    =    50;
small = 1e-10 ;
twopi = 2*pi ;
    

    % --------------------  initialize values   -------------------
    znew  = 0.0;
    dtsec = dtseco;

    if ( abs( dtseco ) > small )
        magro = mag( ro );
        magvo = mag( vo );
        rdotv= dot( ro,vo );

        % -------------  find sme, alpha, and a  ------------------
        sme= ( (magvo^2)*0.5  ) - ( mu /magro );
        alpha= -sme*2.0/mu;

        if ( abs( sme ) > small )
            a= -mu / ( 2.0 *sme );
        else
            a= infinite;
        end
        if ( abs( alpha ) < small )   % parabola
            alpha= 0.0;
        end


        % ------------   setup initial guess for x  ---------------
        % -----------------  circle and ellipse -------------------
        if ( alpha >= small )
            period= twopi * sqrt( abs(a)^3.0/mu  );
            % ------- next if needed for 2body multi-rev ----------
            if ( abs( dtseco ) > abs( period ) )
                % including the truncation will produce vertical lines that are parallel
                % (plotting chi vs time)
                dtsec = rem( dtseco, period );
            end
            xold = sqrt(mu)*dtsec * alpha;
        else
            % --------------------  parabola  ---------------------
            if ( abs( alpha ) < small )
                h = cross( ro,vo );
                magh = mag(h);
                p= magh*magh/mu;
                s= 0.5  * (halfpi - atan( 3.0 *sqrt( mu / (p*p*p) )* dtsec ) );
                w= atan( tan( s )^(1.0 /3.0 ) );
                xold = sqrt(p) * ( 2.0 *cot(2.0 *w) );
                alpha= 0.0;
            else
                % ------------------  hyperbola  ------------------
                temp= -2.0 * mu * dtsec / ...
                    ( a*( rdotv + sign(dtsec)*sqrt(-mu*a)* ...
                    (1.0 -magro*alpha) ) );
                xold= sign(dtsec) * sqrt(-a) *log(temp);
            end
        end

        ktr= 1;
        dtnew = -10.0;
        % conv for dtsec to x units
        tmp = 1.0 / sqrt(mu);

        while ((abs(dtnew*tmp - dtsec) >= small) && (ktr < numiter))
            xoldsqrd = xold*xold;
            znew     = xoldsqrd * alpha;

            % ------------- find c2 and c3 functions --------------
            [c2new, c3new] = findc2c3( znew );

            % ------- use a newton iteration for new values -------
            rval = xoldsqrd*c2new + rdotv*tmp *xold*(1.0 -znew*c3new) + ...
                magro*( 1.0  - znew*c2new );
            dtnew= xoldsqrd*xold*c3new + rdotv*tmp*xoldsqrd*c2new + ...
                magro*xold*( 1.0  - znew*c3new );

            % ------------- calculate new value for x -------------
            temp1 = ( dtsec*sqrt(mu) - dtnew ) / rval;
            xnew = xold + temp1;

            % ----- check if the univ param goes negative. if so, use bissection
            if xnew < 0.0
                xnew = xold*0.5;
            end
  
            ktr = ktr + 1;
            xold = xnew;
        end

        if ( ktr >= numiter )
            errork= 'knotconv';
            fprintf(1,'kep not conv in %2i iter %11.3f \n',numiter, dtseco );
            for i= 1 : 3
                v(i)= 0.0;
                r(i)= v(i);
            end
        else
            % --- find position and velocity vectors at new time --
            xnewsqrd = xnew*xnew;
            f = 1.0  - ( xnewsqrd*c2new / magro );
            g = dtsec - xnewsqrd*xnew*c3new/sqrt(mu);

            for i= 1 : 3
                r(i)= f*ro(i) + g*vo(i);
            end
            magr = mag( r );
            gdot = 1.0  - ( xnewsqrd*c2new / magr );
            fdot = ( sqrt(mu)*xnew / ( magro*magr ) ) * ( znew*c3new-1.0  );
            for i= 1 : 3
                v(i)= fdot*ro(i) + gdot*vo(i);
            end

            temp= f*gdot - fdot*g;
            if ( abs(temp-1.0 ) > 0.00001  )
                errork= 'fandg';
            end

        end  % if
    else
        % ----------- set vectors to incoming since 0 time --------
        for i=1:3
            r(i)= ro(i);
            v(i)= vo(i);
        end
    end





