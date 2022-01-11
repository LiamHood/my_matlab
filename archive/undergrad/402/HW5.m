clear ; close all ; clc ;
%% 1
fprintf( 'Problem 1 \n\n' )
p1()
% Hydrolox
% Find Isp, mass flow, cstar, cf, thrust, nozzle size, storage tank volume

%% 2
fprintf( 'Problem 2 \n\n' )
% Find Isp, At, Ae, mass flow, It, Ab, wp, w0, ws
p2()

%% 3
fprintf( '\n\n Problem 3 \n\n' )
p3()
%% Functions

function p1()

cstarF = @(k,R,Tc) sqrt( k*R*Tc )/( k*sqrt( (2/(k+1))^((k+1)/(k-1)) ) ) ;
CfF = @( k , pc , pe , pa , eps ) sqrt( ((2*k^2)/(k-1)) .* (2/(k+1))^((k+1)/(k-1)) .* (1-(pe./pc).^((k-1)/k)) ) + (( pe - pa )./pc ).*eps ;
% Find Isp, mass flow, cstar, cf, thrust, nozzle size, storage tank volume
fprintf( 'I used the included graphs to find the initial values for combustion \n' )
fprintf( 'temperature, the gamma, the molecular weight, and the Isp. \n' )
names = [ "Stoichiometric Ratio \n" , "Max Isp \n" , "Balance of the Two \n" ] ;

    dv  =   500     ; % m/s
    m   =   20000   ; % kg
    tb  =   10*60   ; % s
    pc  =   2e6     ; % Pa
    mrv  =   [ 8 , 3.8 , 6.2 ] ; % Mixture ratio
    Tcv  =   [ 3500 , 2800 , 3450 ] ; % combustion temp
    kv   =   [ 1.2 , 1.22 , 1.205 ] ;
    Mv   =   [ 16 , 9.5 , 13.7 ]*1e-3 ;
    Rbar =   8.31446261815324 ;
    eps = 20 ;
    g0 = 9.81 ;
    
    figure
    for ii = 1:3 
        mr = mrv(ii) ;
        Tc = Tcv(ii) ;
        k = kv(ii) ;
        M = Mv(ii) ;
        n = 1e3 ;

        cstar(ii) = cstarF( k , Rbar/M , Tc ) ;
        peFull = linspace( 1e3 , 1e4 , n ) ;
        pcFull = linspace( 1e6 , 1e7 , n ) ;
        for jj = 1:n
            for kk = 1:n
                eps(jj,kk) = 1/( ((k+1)/2)^(1/(k-1)) .* (peFull(jj)./pcFull(kk)).^(1/k) ...
                    .* sqrt( ((k+1)./(k-1)).*(1-(peFull(jj)./pcFull(kk)).^((k-1)/k) ) ) ) ;
            end
        end
        Cf = CfF( k , pcFull , peFull , 0 , eps ) ;
%         surf( pcFull , peFull , eps )
%         xlabel( 'Chamber Pressure' )
%         ylabel( 'Exit Pressure' )
%         zlabel( 'Area Ratio' )
%         hold on
        
        fprintf( 'Lower Chamber pressure gives a low area ratio \n' )
        pc = 1e6 ;
        epsFull = 1./( ((k+1)/2)^(1/(k-1)) .* (peFull./pc).^(1/k) ...
                    .* sqrt( ((k+1)./(k-1)).*(1-(peFull./pc).^((k-1)/k) ) ) ) ;
        Cf = CfF( k , pc , peFull , 0 , epsFull ) ;
        plot( epsFull , Cf )
        ylabel( 'Coefficient of Thrust' )
        xlabel( 'Area Ratio' )
        hold on
        
        eps = 30 ;
        pe(ii) = peFull( find( epsFull >= eps & epsFull < ( eps + eps*2e-3 ) ) ) ;
        Cfreal(ii) = Cf( find( epsFull >= eps & epsFull < ( eps + eps*2e-3 ) ) ) ;
        c(ii) = cstar(ii)*Cfreal(ii) ;
        Isp(ii) = c(ii)/g0 ;
        
        
        
%         thrust(ii) = Cf(ii)*pc'*At ;
%         mdot(ii) = thrust(ii)/(Isp(ii)*g0) ;
%         mp(ii) = mdot(ii)*tb ;
        mp(ii) = m*exp( dv/c(ii) ) - m ;
        mdot(ii) = mp(ii)/tb ;
        thrust(ii) = c(ii)*mdot(ii) ;
        At(ii) = thrust(ii)./( Cfreal(ii).*pc ) ;
        dt(ii) = sqrt( At(ii)/pi ) ;
        Ae(ii) = At(ii)*eps ;
        de(ii) = sqrt( Ae(ii)/pi ) ;
        vol(ii) = ( mp(ii) / ( .071 * ( 1/1000 ) * ( 1000/1 ) ) ) ;
        
        
        
fprintf( names(ii) )
fprintf( 'The Specific Impulse is: %f s \n' , Isp(ii) )
fprintf( 'The Coefficient of Thrust is: %f \n' , Cfreal(ii) )
fprintf( 'The mass of propellant is: %f kg \n' , mp(ii) )
fprintf( 'The mass flow rate is: %f kg/s \n' , mdot(ii) )
fprintf( 'The throat diameter is: %f m \n' , dt(ii) )
fprintf( 'The nozzle exit diameter is: %f m \n' , de(ii) )
fprintf( 'The thrust is: %f N \n' , thrust(ii) )
fprintf( 'The fuel volume is: %f L \n \n' , vol(ii) )


    end
  fprintf( 'Fuel density from http://www.astronautix.com/l/loxlh2.html \n \n' )      
end

function p2()
cstarF = @(k,R,Tc) sqrt( k*R*Tc )/( k*sqrt( (2/(k+1))^((k+1)/(k-1)) ) ) ;
CfF = @( k , pc , pe , pa , eps ) sqrt( ((2*k^2)/(k-1)) * (2/(k+1))^((k+1)/(k-1)) * (1-(pe/pc)^((k-1)/k)) ) + (( pe - pa )/pc )*eps ;
CfsF = @( k , pc , pa ) sqrt( ((2*k^2)/(k-1)) .* (2/(k+1))^((k+1)/(k-1)) .* (1-(pa./pc).^((k-1)/k)) ) ;

    p0      =   101325 ; % Pa
    eff     =   .98     ;
    mpCorr  =   1.04    ;
    It2w    =   143     ; % s
    gamma   =   1.26    ;
    Tcbad      =   2700    ; % Fahrenheit
    Tc         =   1755.372 ; % Kelvin
    rdotbbad   =   0.10    ; % in/s @1000 psi
    rdotb      =   rdotbbad*2.54 ; % cm/s
    cstarbad   =   4000    ; % ft/s
    cstar      =   cstarbad*.3048 ; % m/s
    rhopbad    =   .056    ; % lb/in^3
    rhop       =   rhopbad*0.453592*(1/.0254)^3 ; % kg/m^3
    Mbad       =   22      ; % lbm/(lb*mole)
    M          =   Mbad*(0.453592/453.59) ;
    Ftbad      =   2000    ; % lbf
    Ft         =   4.44822*Ftbad ; % N
    tb      =   10      ; % s
    Pcbad     =   1000    ; % psia
    pc          =  Pcbad*6894.76 ; % Pa
    Tobad      =   70      ; % Fahrenheit
    To         = 294.261 ; % kelvin
    prop    =   "Ammonium Nitrate" ;
    
    Cf = CfsF(gamma,pc,p0) ;
    
    Isp = cstar*Cf*eff/9.81 ;
    At = Ft./( Cf*eff.*pc ) ;
    Ae = At./( ((gamma+1)/2)^(1/(gamma-1)) .* (p0./pc).^(1/gamma) .* sqrt( ((gamma+1)/(gamma-1)).*(1-(p0./pc).^((gamma-1)/gamma) ) ) ) ;
    mdot = pc.*At./cstar ;
    It = tb*Ft ;
    Ab = mdot/(rhop*rdotb) ;
    mprop = (mdot.*tb)*1.04 ;
    wprop = mprop*9.81 ;
    w = It/It2w ;
    
    fprintf( 'The Solid Rocket Motor characteristics are: \n' )
    fprintf( 'Specific Impulse: %f s \n' , Isp )
    fprintf( 'Throat Area: %f m^2 \n' , At )
    fprintf( 'Nozzle Exit Area: %f m^2 \n' , Ae )
    fprintf( 'Mass Flow Rate: %f kg/s \n' , mdot )
    fprintf( 'Total Impulse: %f N*s \n' , It )
    fprintf( 'Burning Area: %f m^2 \n' , Ab )
    fprintf( 'Propellant Weight: %f N \n' , wprop )
    fprintf( 'Gross Weight: %f N \n' , w )
    

    
end

function p3()
cstarF = @(k,R,Tc) sqrt( k*R*Tc )/( k*sqrt( (2/(k+1))^((k+1)/(k-1)) ) ) ;
CfsF = @( k , pc , pa ) sqrt( ((2*k^2)/(k-1)) .* (2/(k+1))^((k+1)/(k-1)) .* (1-(pa./pc).^((k-1)/k)) ) ;

    Ft  =   7500    ; % N
    tb  =   100     ; % s
    alt =   10      ; % km
    p0  =   26.5e3  ; % Pa
    g0  =   9.81    ; % m/s^2
    
    %D
    rhop = 1680 ;
    gamma = 1.2 ;
    cstar = 1511 ;
    Tc = 3288 ;
    pc = linspace( 2.7 , 10 , 1e3 )*1e6 ;
    n = .39 ;
    a = .41 ;
    
    
    Cf = CfsF( gamma , pc , p0 ) ;
    At = Ft./( Cf.*pc ) ;
    dt = 2*sqrt( At./pi ) ;
    figure
    plot( pc*1e-6 , dt*1e2 )
    xlabel( 'Chamber Pressure [MPa]' )
    ylabel( 'Throat Diameter [cm]' )
    
    Ae = At./( ((gamma+1)/2)^(1/(gamma-1)) .* (p0./pc).^(1/gamma) .* sqrt( ((gamma+1)/(gamma-1)).*(1-(p0./pc).^((gamma-1)/gamma) ) ) ) ;
    de = 2*sqrt( Ae./pi ) ;
    figure
    plot( pc*1e-6 , de*1e2 )
    xlabel( 'Chamber Pressure [MPa]' )
    ylabel( 'Nozzle Exit Diameter [cm]' )
    
    rb = a.*(pc*1e-6).^n ;
%     figure
%     plot( pc*1e-6 , rb )
%     xlabel( 'Chamber Pressure [MPa]' )
%     ylabel( 'Fuel Grain Burn Rate [cm/s]' )
    
    mdot = pc.*At./cstar ;
    mprop = mdot.*tb ;
    figure
    plot( pc*1e-6 , mprop )
    xlabel( 'Chamber Pressure [MPa]' )
    ylabel( 'Mass of Propellant [kg/s]' )
    
    rgrain = rb.*tb ;
    figure
    plot( pc*1e-6 , rgrain )
    xlabel( 'Chamber Pressure [MPa]' )
    ylabel( 'Fuel Grain Thickness [cm]' )
    
    Isp = Ft./(mdot.*g0) ;
    figure
    plot( pc*1e-6 , Isp )
    xlabel( 'Chamber Pressure [MPa]' )
    ylabel( 'Specific Impulse [s]' )
    
    fprintf( 'The solid rocket motor should be designed around the maximum \n' )
    fprintf( 'chamber pressure that the fuel grain is designed for. This is \n' )
    fprintf( 'a pressure of %f MPa. This will need some extra mass for containing \n' , pc(end)*1e-6 )
    fprintf( 'the pressure but it save mass with a smaller nozzle and propellant. \n' )
    fprintf( 'For this chamber pressure the characteristics are as follow: \n' )
    fprintf( 'Chamber Pressure: %f MPa \n' , pc(end)*1e-6 )
    fprintf( 'Throat Diameter: %f cm \n' , dt(end) )
    fprintf( 'Nozzle Exit Diameter: %f cm \n' , de(end) )
    fprintf( 'Grain Thickness: %f cm \n' , rgrain(end) )
    fprintf( 'Mass of Propellant: %f kg \n' , mprop(end) )
    fprintf( 'Specific Impulse: %f s \n' , Isp(end) )
end


