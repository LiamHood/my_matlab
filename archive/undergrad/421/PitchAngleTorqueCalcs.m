clear ; close all ; clc ;
% Set the values for the configurations with different beam lengths
L = [ 1 , 2 , 5 , 7 , 10 , 15 , 20 ] ;
zoff = [ .2142 , .2807 , .5400 , .7628 , 1.1719 , 2.0533 , 3.1842 ] ;
COM = [ zeros( 2 , length( zoff ) ) ; -zoff ] ;
Iz = [ .1271 , .1259 , .1221 , .1197 , .1159 , .1097 , .1034 ] ;
Ix = [ 1.5571 , 5.1660 , 32.818 , 68.2014 , 152.0170 , 388.8188 , 766.7421 ] ;
Iy = [ 1.5571 , 5.1660 , 32.818 , 68.2015 , 152.0172 , 388.8191 , 766.7426 ] ;

% Set orbital parameters. Without the aerospace blockset and my own state
% to coes functions I don't want to use actual state vectors so I am using
% vectors of the same magnitude as the ISS but in a different plane. 

r = [ 0 ; 6.7927e+06 ; 0 ] ;
v = [ 7.6567e+03 ; 0 ; 0 ] ;
ang = linspace( 0 , 90 , 1e3 ) ;

% Calculate gravity gradient over a variety of angles and all lengths
    % the rotation about the y axis brings the beam to point in the nadir
    % direction
for jj = 1:length(ang)
    Cbg = roty( ang(jj) )*rotx( 90 ) ;
    for ii = 1:length(L)
        Tgg(:,ii,jj) = GravGradFun( Cbg , Ix(ii) , Iy(ii) , Iz(ii) , r ) ;
        Tggmag(ii,jj) = norm( Tgg(:,ii,jj) ) ;
    end
end

figure 
semilogy( ang , Tggmag(1,:) , ang , Tggmag(2,:) , ang , Tggmag(3,:) , ang , Tggmag(4,:) , ang , Tggmag(5,:) , ang , Tggmag(6,:) , ang , Tggmag(7,:) )
legend( '1 meter' , '2 meters' , '5 meters' , '7 meters', '10 meters' , '15 meters' , '20 meters' )
title( 'Gravity Gradient Torque Magnitude' )
ylabel( 'Torque (Nm)' )
xlabel( 'Pitch Angle (degrees)' )


% Calculate atmospheric drag torque over a variety of angles and all lengths
    % the rotation about the y axis brings the beam to point in the nadir
    % direction
for jj = 1:length(ang)
    Cbg = roty( ang( jj ) )*rotx( 90 ) ;
    for ii = 1:length(L)
        % areas of all faces
        abzp = .04 ;
        abzn = .04 - .0234*.04684 ;
        abyp = .06 ;
        abyn = .06 ;
        abxp = .06 ;
        abxn = .06 ;
        alyp = .04684*L(ii) ;
        alyn = .04684*L(ii) ;
        alxp = .0234*L(ii) ;
        alxn = .0234*L(ii) ;
        aezp = .01 - .0234*.04684 ;
        aezn = .01 ;
        Areas = [ abzp , abzn , abyp , abyn , abxp , abxn , alyp , alyn , alxp , alxn , aezp , aezn ] ;
        % positions of all faces
        pbzp = [ 0 ; 0 ; 0 ] ;
        pbzn = [ 0 ; 0 ; -.3 ] ; 
        pbyp = [ 0 ; .1 ; -.15 ] ;
        pbyn = [ 0 ; -.1 ; -.15 ] ;
        pbxp = [ .1 ; 0 ; -.15 ] ;
        pbxn = [ -.1 ; 0 ; -.15 ] ;
        plyp = [ 0 ; .04684/2 ; -.15-L(ii)/2 ] ;
        plyn = [ 0 ; -.04684/2 ; -.15-L(ii)/2 ] ;
        plxp = [ .0234/2 ; 0 ; -.15-L(ii)/2 ] ;
        plxn = [ .0234/2 ; 0 ; -.15-L(ii)/2 ] ;
        pezp = [ 0 ; 0 ; -L(ii) ] ;
        pezn = [ 0 ; 0 ; -L(ii) ] ;
        Positions = [ pbzp , pbzn , pbyp , pbyn , pbxp , pbxn , plyp , plyn , plxp , plxn , pezp , pezn ] - COM(:,ii) ;
        % normal vectors of faces
        nx = [ 1 ; 0 ; 0 ] ;
        ny = [ 0 ; 1 ; 0 ] ;
        nz = [ 0 ; 0 ; 1 ] ;
        Norms = [ nz , -nz , ny , -ny , nx , -nx , ny , -ny , nx , -nx , nz , -nz ] ;
        %Find drag Torque
        Tatmo(:,ii,jj) = DragFun( Cbg , v , r , Areas , Positions , Norms ) ;
        Tatmomag(ii,jj) = norm( Tatmo(:,ii,jj) ) ;
        % Find Solar radiation pressure
        s = [ 1 ; 0 ; 0 ] ;
        Tsrp(:,ii,jj) = SolarPressFun( Cbg , r , s , Areas , Positions , Norms ) ;
        Tsrpmag(ii,jj) = norm( Tsrp(:,ii,jj) ) ;
        
        Torque = Tsrp + Tatmo + Tsrp ;
        Torquemag( ii , jj ) = norm( Torque(:,ii,jj) ) ; 
    end
end

figure 
semilogy( ang , Tatmomag(1,:) , ang , Tatmomag(2,:) , ang , Tatmomag(3,:) , ang , Tatmomag(4,:) , ang , Tatmomag(5,:) , ang , Tatmomag(6,:) , ang , Tatmomag(7,:) )
legend( '1 meter' , '2 meters' , '5 meters' , '7 meters', '10 meters' , '15 meters' , '20 meters' )
title( 'Atmospheric Drag Torque Magnitude' )
ylabel( 'Torque (Nm)' )
xlabel( 'Pitch Angle (degrees)' )

figure 
semilogy( ang , Tsrpmag(1,:) , ang , Tsrpmag(2,:) , ang , Tsrpmag(3,:) , ang , Tsrpmag(4,:) , ang , Tsrpmag(5,:) , ang , Tsrpmag(6,:) , ang , Tsrpmag(7,:) )
legend( '1 meter' , '2 meters' , '5 meters' , '7 meters', '10 meters' , '15 meters' , '20 meters' )
title( 'Solar Radiation Pressure Torque Magnitude' )
ylabel( 'Torque (Nm)' )
xlabel( 'Pitch Angle (degrees)' )

% Gravity Gradient Torque is trying to keep the beam pointed nadir while
% the atmospheric drag torque is trying to keep the beam in the wake
% direction at some angle they will balance

% find balance
for ii = 1:length(L)
    diff(:,ii,:) = Tgg(:,ii,:) + Tatmo(:,ii,:) + Tsrp(:,ii,:)  ;
    sdiff = size( diff ) ;
        diffreal(ii,:) = diff( 2 , ii , : ) ;
        for jj = 1:length( ang )
            diffmag( ii,jj )= norm( diffreal( ii , jj ) ) ;
        end
    [ mintorque , indmin ] = min( diffmag( ii , 1:end/2 ) ) ;
    balanceAngle( ii ) = ang( indmin ) ;
end

figure
semilogy( ang , diffmag(1,:) , ang , diffmag(2,:) , ang , diffmag(3,:) , ang , diffmag(5,:) , ang , diffmag(end,:) )
legend( '1 meter' , '2 meters' , '5 meters' ,'10 meters' , '20 meters' , 'Location' , 'south' )
title( 'Torque Magnitude' )
ylabel( 'Torque (Nm)' )
xlabel( 'Pitch Angle (degrees)' )

figure 
plot( L , balanceAngle )
title( 'Pitch Angle that Torque balances at' )
ylabel( 'Pitch Angle (degrees)' )
xlabel( 'Length of Beam (m)' )







