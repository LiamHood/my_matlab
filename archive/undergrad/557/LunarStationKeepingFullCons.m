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

load( 'YoD_FC.mat' ) ;
raanDriftAvg = zeros( 1 , 38 ) ;
for ii = 2:38
    raanDriftAvg(ii) = mean( raanP(ii,:) - raanP(1,:) ) ;
end

for plane = 1:4
    for d = 2:38
    % set arrival / nominal orbit
%     incNom = 55 ;
    incNom = incP(d,plane) ;
    aopNom = 0 ; 
    eccNom = 1e-9 ;
    hNom = hP(1,plane) ;
    thetaS = 0 ;
    thetaEnd = 360 ;
    thetaNom = linspace( thetaS, thetaEnd , m ) ;
    % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
    rNom = zeros( 3 , m ) ;
    vNom = zeros( 3 , m ) ;
%     raanNom = raanP(1,plane) + raanDriftAvg(d) ;
    raanNom = raanP(d,plane) ;
        for ii = 1:m
            [ rNom(:,ii) , vNom(:,ii) ] = coes2state( [ hNom , incNom*d2r , eccNom , raanNom*d2r , aopNom(1)*d2r , thetaNom(ii)*d2r ] , mu ) ;
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

        [ min_perCol , rows ] = min( deltav ) ;
        [ mindv(d) , col ] = min( min_perCol ) ;
        row = rows( col ) ;

        tol = 1e-8 ;
        [ t , rv ] = TwoBody( [ 0 , 1.1*T ] , rNom(:,col) , vNom(:,col) , mu , tol ) ;
        [ tPer , rvPer ] = TwoBody( [ 0 , 1.1*T ] , rPer , vPer , mu , tol ) ;
        [ tTran , rvTran ] = TwoBody( [ 0 , tf( row , col ) ] , rPer , va(:,row,col) , mu , tol ) ;

        dir = [ "x" ; "y" ; "z" ] ;
        burn = va(:,row,col) - vPer ;
        burn2 = vNom(:,col) - vb(:,row,col) ;
        dispva = [ burn(1) , dir(1) , burn(2) , dir(2) , burn(3) , dir(3) ] ;
        dispvb = [ burn2(1) , dir(1) , burn2(2) , dir(2) , burn2(3) , dir(3) ] ;
        fprintf( 'The delta-v is %f m/s \n' , deltav(row,col)*1e3 )
        fprintf( 'The time of transfer is %f hours \n' , tf(row,col)*(1/(60*60)) )
        daysnum(d) = (d-1)*10 ;
        dvyear(plane,d) = mindv(d)*(365/daysnum(d)) ;
        thetaEndPP(plane) = thetaNom( col ) ;
    end
end

for plane = 5:32
    for d = 2:38
    trueplane = plane - 4*floor( plane/4 ) ;
    if trueplane == 0
        trueplane = 4 ;
    end
    sinplane = floor( (plane-1)/4) + 1 ;
    % set arrival / nominal orbit
%     incNom = 55 ;
    incNom = incP( d,plane );
    aopNom = 0 ; 
    eccNom = 1e-9 ;
    hNom = hP(1,plane) ;
    thetaNom = thetaEndPP( trueplane ) + ( pi/4 )*sinplane ;
    % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
    rNom = zeros( 3 , m ) ;
    vNom = zeros( 3 , m ) ;
%     raanNom = raanP(1,plane) + raanDriftAvg(d) ;
    raanNom = raanP( d,plane ) ;
    clear rNom ;
    clear vNom ;
    [ rNom(:,1) , vNom(:,1) ] = coes2state( [ hNom , incNom(1)*d2r , eccNom , raanNom*d2r , aopNom*d2r , thetaNom*d2r ] , mu ) ;

        % set beginning / perturbed orbit
    %     incPer = [55.2818133233726;56.3257284620972;56.9635291085618;56.8906611133049;57.3164432916295;58.3723838156408;58.9685428999927;59.0041674165460;59.3826542638042;60.434174268500510] ;
        
        smaPer = (hP(d,plane)^2)/(mu*(1-eccP(d,plane)^2)) ;
        T = 2*pi*sqrt( smaPer^3 / mu ) ;
        % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
        [ rPer , vPer ] = coes2state( [ hP(d,plane) , incP(d,plane)*d2r , eccP(d,plane) , raanP(d,plane)*d2r , aopP(d,plane)*d2r , thetaP(d,plane)*d2r ] , mu ) ;
        
        
        % find min time
        clear tp ;
        clear tmin ;
        [ ~ , ~ , ~ , tmin(:) , tp ] = Lambert_MinEnergy( rPer , rNom , mu ) ;
        

        % set time of flight
        clear tf
        tf(:) = linspace( tp , T*2 , n ) ;

        % use lamberts solver
        clear deltav ;
        va = zeros( 3 , n ) ;
        vb = zeros( 3 , n ) ;
        deltav = zeros( 1,n ) ;
        collision = zeros( 1,n ) ;
        for ii = 1:n
                [ va(:,ii) , vb(:,ii) , deltav(ii) , collision(ii) ] = Targeting_FF( rPer , rNom , vPer , vNom , tf(ii) , mu , Rcb ) ;
                if collision(ii) == 1 
                    deltav(ii) = 10 ;
                end
        end

        [ mindv(d) , rows ] = min( deltav ) ;

        tol = 1e-8 ;
        [ t , rv ] = TwoBody( [ 0 , 1.1*T ] , rNom , vNom , mu , tol ) ;
        [ tPer , rvPer ] = TwoBody( [ 0 , 1.1*T ] , rPer , vPer , mu , tol ) ;
        [ tTran , rvTran ] = TwoBody( [ 0 , tf( rows ) ] , rPer , va(:,rows) , mu , tol ) ;

        dir = [ "x" ; "y" ; "z" ] ;
        burn = va(:,rows) - vPer ;
        burn2 = vNom(:) - vb(:,rows) ;
        dispva = [ burn(1) , dir(1) , burn(2) , dir(2) , burn(3) , dir(3) ] ;
        dispvb = [ burn2(1) , dir(1) , burn2(2) , dir(2) , burn2(3) , dir(3) ] ;
        fprintf( 'The delta-v is %f m/s \n' , deltav(rows)*1e3 )
        fprintf( 'The time of transfer is %f hours \n' , tf(rows)*(1/(60*60)) )
        
        daysnum(d) = (d-1)*10 ;
        dvyear(plane,d) = mindv(d)*(365/daysnum(d)) ;
    end
end


    
    plot( daysnum , dvyear )
    axis([0,370,0,1])
    hold on


legend( 'Plane 1' , 'Plane 2' , 'Plane 3' , 'Plane 4' )