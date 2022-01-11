%% Final Project
clear ; close all ; clc ;
mu = 398600 ;
tol = 1e-10 ;
forces = 'gravity , J2 , J3 , drag , srp , nbody ' ;

%% Flying laptop
% h , inc , ecc , RAAN , omega , theta , a , rp , ra 
m = 120 ;
areaD = .7*.9 ;
areaS = 1 ;
A = areaS ;
d2r = pi/180 ;
epoch = [ 2019 , 12 , 1 ] ;
JDo = juliandate( epoch ) ;
tspan = [ 0 , 120*24*60*60 ] ;
state0 = 1e3*[  -6.064781955136040 ; -2.105667897893976;  -2.698157301527066;  -0.003098601440846;   0.000065711340527;   0.006906944507578 ] ;
sf = [ 1 , 1 , 1 ] ;
dt = 100 ;
COESo = state2COE( state0 , mu ) ;
Po = 2*pi.*sqrt( COESo(7).^3 ./ mu ) ;
tspan2 = [ 0 , 120*24*60*60 + Po/2 ] ;
% [ tp , COESp ] = VoP( tspan , state0 , mu , tol , forces , A , m , epoch , JDo ) ;
[ tp , tpc , rvp , COESpc] = Cowell( tspan , state0 , mu , tol , forces , A , m , sf , epoch , JDo ) ;
[ t , tadj , rv , COES ] = TwoBody( tspan , state0 , mu , tol , sf ) ;
        Ppc = 2*pi.*sqrt( COESpc(:,7).^3 ./ mu )*(1/3600) ;
%         Pp = 2*pi.*sqrt( COESp(:,7).^3 ./ mu )*(1/3600) ;
        P = 2*pi.*sqrt( COES(:,7).^3 ./ mu )*(1/3600) ;
        s2d = 1/(24*60*60) ;
        r2d = 180/pi ;
        figure
        subplot( 2 , 1 , 1 )
        hold on
        plot( tadj*s2d , COES(:,2)*r2d )
%         plot( tp*s2d , COESp(:,2)*r2d )
        plot( tpc*s2d , COESpc(:,2)*r2d )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'Inclination [degrees]' )
        legend( 'Unperturbed' , 'Perturbed' )

        subplot( 2 , 1 , 2 )
        hold on
        plot( tadj*s2d , COES(:,3) )
%         plot( tp*s2d , COESp(:,3) )
        plot( tpc*s2d , COESpc(:,3) )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'Eccentricity' )
        legend( 'Unperturbed' , 'Perturbed' )

        figure
        subplot( 2 , 1 , 1 )
        hold on 
        plot( tadj*s2d , COES(:,4)*r2d )
%         plot( tp*s2d , COESp(:,4)*r2d )
        plot( tpc*s2d , COESpc(:,4)*r2d )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'RAAN [degrees]' )
        legend( 'Unperturbed' , 'Perturbed' )

        subplot( 2 , 1 , 2 )
        hold on
        plot( tadj*s2d , COES(:,5)*r2d )
%         plot( tp*s2d , COESp(:,5)*r2d )
        plot( tpc*s2d , COESpc(:,5)*r2d )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'Argument of Perigee [degrees]' )
        legend( 'Unperturbed' , 'Perturbed' )

        figure
        subplot( 2 , 1 , 1 )
        hold on
        plot( tadj*s2d , COES(:,7) )
%         plot( tp*s2d , COESp(:,7) )
        plot( tpc*s2d , COESpc(:,7) )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'Semi Major Axis [km]' )
        legend( 'Unperturbed' , 'Perturbed' )
        
        subplot( 2 , 1 , 2 )
        hold on
        plot( tadj*s2d , P )
%         plot( tp*s2d , Pp )
        plot( tpc*s2d , Ppc )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'Period []' )
        legend( 'Unperturbed' , 'Perturbed' )

        figure
        subplot( 2 , 1 , 1 )
        hold on
        plot( tadj*s2d , COES(:,8) )
%         plot( tp*s2d , COESp(:,7) )
        plot( tpc*s2d , COESpc(:,8) )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'Rp [km]' )
        legend( 'Unperturbed' , 'Perturbed' )
        
        subplot( 2 , 1 , 2 )
        hold on
        plot( tadj*s2d , COES(:,8) )
%         plot( tp*s2d , Pp )
        plot( tpc*s2d , COESpc(:,9)  )
        hold off
        xlabel( 'Time [days]' )
        ylabel( 'Ra [km]' )
        legend( 'Unperturbed' , 'Perturbed' )
% 
% [ v1 , v2 ] = LambertsBook( rvp(end,1:3) , rv(end,1:3) , Po/2 , mu , tol , 1 ) ;
% [ tburn , tadjburn , rvburn , COESburn ] = TwoBody( [0 , Po/2] , state0 , mu , tol , sf ) ;
% figure
% hold on
% plot3( rvp(:,1:3) )
% plot3( rv(:,1:3) )
% plot3( rvburn(:,1:3) )

    

% %% Viasat - 1
% % h , inc , ecc , RAAN , omega , theta , a , rp , ra 
% mass = 6740 ;
% areaD = 3.409*8.733 ;
% areaS = 1 ;
% d2r = pi/180 ;
% 
% inc0 = 0.0186*d2r ;
% ecc0 = .000226 ;
% raan0 = 265.0753*d2r ;
% omega0 = 311.649*d2r ;
% 
% 
% COES0 = [ h0 , inc0 , ecc0 , raan0 , omega0 , ta0 ] ;
% [ r0 , v0 ] = coes2state( COES0 , mu ) ;
% [ t , COES ] = VoP( tspan , [ r0 ; v0 ] , mu , tol , forces , A , m , epoch ) ;
