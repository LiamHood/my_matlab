%% Homework 6 
% Aero 402
% Liam Hood
function HW6()
clear ; close all ; clc ;
pt = 'Problem number %u \n' ;

%% 1
fprintf( pt , 1 ) 
HW6_P1()
fprintf( ' \n' )

%% 2 
fprintf( pt , 2 ) 
HW6_P2()
fprintf( ' \n' )

%% 3 
fprintf( pt , 3 ) 
HW6_P3()
fprintf( ' \n' )

%% Functions
    function HW6_P1()
       
        dg = 20e-2 ;
        Lo = 50e-2 ;
        dp = 4e-2 ;
        a = 0.0018 ;
        n = 0.3 ;
        rhoFUEL = 920 ;
        mox = 4 ;
        tspan = [ 0 , 8 ] ;
        
        
        
        opts = odeset( 'AbsTol' , 1e-8 , 'RelTol' , 1e-8 , 'Events' , @BurnOut ) ;
        [ tf , rf ] = ode45( @RDot , tspan , dp/2 , opts , a , mox , n , dg/2 ) ;
        figure
        plot( tf , rf*2*1e2 )
        xlabel( 'Time [s]' )
        ylabel( 'Port Diameter [cm]' )
        
        rrate = RDot( tf , rf , a , mox , n ) ;
        figure
        plot( tf , rrate*1e2 )
        xlabel( 'Time [s]' )
        ylabel( 'Regression Rate [cm/s]' )
        
        mdotf = 2*pi*rhoFUEL.*rf.*Lo.*rrate ;
        mdot = mdotf + mox ;
        
        figure
        axis( [ 0 2 0 6 ] )
        plot( tf , mox*ones(length(tf),1) , tf , mdotf , tf , mdot ) 
        xlabel( 'Time [s]' )
        ylabel( 'Mass Flow Rate [kg/s]' )
        legend( 'Oxidizer' , 'Fuel' , 'Total' )
    end

    function HW6_P2()
        CfF = @( k , pc , pe , pa , eps ) sqrt( ((2*k^2)/(k-1)) * (2/(k+1))^((k+1)/(k-1)) .* (1-(pe./pc).^((k-1)/k)) ) + (( pe - pa )./pc ).*eps ;
        
        Fi = 12e6 ;
        pci = 4.8e6 ;
        o2f = 2 ;
        k = 1.52 ;
        cstar = 1800 ;
        dg = 3.8 ;
        Lo = 40 ;
        dp = .59 ;
        db = .34 ;
        a = 0.0018 ;
        n = 0.07 ;
        rhoFUEL = 920 ;
        rhoOX = 1170 ;
        mox = 2844 ;
        tspan = [ 0 , 120 ] ;
        ports = 7 ;
        At = 1.6 ;
        pa = 80 ;
        eps = 6.67 ;
        
        opts = odeset( 'AbsTol' , 1e-8 , 'RelTol' , 1e-8 , 'Events' , @BurnOut ) ;
        [ tf , rf ] = ode45( @RDot , tspan , dp/2 , opts , a , mox , n , (db*2+dp)/2 ) ;
        figure
        plot( tf , rf*2 )
        xlabel( 'Time [s]' )
        ylabel( 'Port Diameter [m]' )
        
        rrate = RDot( tf , rf , a , mox , n ) ;
        figure
        plot( tf , rrate*1e2 )
        xlabel( 'Time [s]' )
        ylabel( 'Regression Rate [cm/s]' )
        fprintf( 'Rocket actually burns out all of the fuel before 120s \n' )
        fprintf( 'It burns out at %f seconds \n' , tf(end) )
        
        mdotf = 2*pi*rhoFUEL*rf(end)*Lo*rrate(end) ;
        mdot = mdotf+mox  ;
        fprintf( 'The O/F ratio is %f \n' , mox/mdotf )
        fprintf( 'The mass flow rate of propellant is %f [kg/s] \n' , mdot )
        
        
        pc = mdot*cstar/At ;
        fprintf( 'The chamber pressure is %f [MPa] \n' , pc*1e-6 )
        
        nn = 1e4 ;
        peFull = linspace( 0 , 1e5 , nn ) ;
        epsfull = 1./( ((k+1)/2)^(1/(k-1)) .* (peFull./pc).^(1/k) ...
                    .* sqrt( ((k+1)./(k-1)).*(1-(peFull./pc).^((k-1)/k) ) ) ) ;
        pe = peFull( find( epsfull >= 6.67 & epsfull <= 6.671 ) ) ;
        
        Cf = CfF( k , pc , pe , pa , eps ) ;
        Ff = Cf*pc*At ;
        fprintf( 'The thrust at 50 km is %f [MN] \n' , Ff*1e-6 )
        
        Isp = Ff/(mdot*9.81) ;
        fprintf( 'The Isp at 50 km is %f [s] \n' , Isp )
    end

    function HW6_P3()
        % find thrust , Isp , It
        AreaRatio = @(k,peFull,pc)1./( ((k+1)/2)^(1/(k-1)) .* (peFull./pc).^(1/k).* sqrt( ((k+1)./(k-1)).*(1-(peFull./pc).^((k-1)/k) ) ) ) ;
        CfF = @( k , pc , pe , pa , eps ) sqrt( ((2*k^2)/(k-1)) * (2/(k+1))^((k+1)/(k-1)) .* (1-(pe./pc).^((k-1)/k)) ) + (( pe - pa )./pc ).*eps ;
        cstarF = @(k,R,Tc) sqrt( k*R*Tc )/( k*sqrt( (2/(k+1))^((k+1)/(k-1)) ) ) ;

        k = 1.2 ;
        R = 259 ;
        Tc = 3000 ;
        pc = 200e3 ;
        At = .5*1e-4 ;
        eps = 15 ;
        mp = 10 ;
        
        nn = 1e4 ;
        peFull = linspace( 1e3 , 1e4 , nn ) ;
        epsfull = 1./( ((k+1)/2)^(1/(k-1)) .* (peFull./pc).^(1/k) ...
                    .* sqrt( ((k+1)./(k-1)).*(1-(peFull./pc).^((k-1)/k) ) ) ) ;
        pe = peFull( find( epsfull >= eps & epsfull <= eps+.005 ) ) ;
        Cf = CfF( k , pc , pe , 0 , eps ) ;
        F = Cf*pc*At ;
        fprintf( 'The thrust is %f [N] \n' , F )
        
        cstar = cstarF(k,R,Tc) ;
%         mdot = At*pc/cstar ;
        Isp = cstar*Cf/9.81 ;
        fprintf( 'The specific impulse is %f [s] \n' , Isp )
        
        It = Isp*mp*9.81 ;
        fprintf( 'The total impulse is %f [kN*s] \n' , It*1e-3 )
        
        
    end

%% Function
    function rdot = RDot(t,r,a,mox,n,rmax) 
        rdot = a.*(mox./(r.^2.*pi)).^n ;
    end
    function [ value , isterminal , direction ] = BurnOut( t , r , a , mox , n , rmax )
        value = r - rmax ;
        isterminal = 1 ;
        direction = 0 ;
    end
end
