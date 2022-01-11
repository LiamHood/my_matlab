%% Homework 7
% Aero 402
% Liam Hood
function HW7()
clear ; close all ; clc ;
pt = 'Problem number %u \n' ;

% 2
fprintf( pt , 2 ) 
HW7_P2()
fprintf( ' \n' )


%% Problems
    function HW7_P2() 
        FepF = @(Ft,alpha,Ib,Mi,Vb,q) Ft.*alpha.*Ib.*sqrt( ( 2.*Mi.*Vb )./q );
        IspepF = @(Ft,alpha,mue,Mi,Vb,q) ( ( Ft.*alpha.*mue )./9.81 )*sqrt( ( 2.*q.*Vb )./Mi ) ;
        
        d = .15 ;
        Mxe = 2.18e-25 ;
        Pin = 3e3 ;
        Vb = 800 ;
        Ib = 2.2 ;
        theta = 10 ;
        mue = .85 ;
        massp = 20 ;
        q = 1.60217662e-19 ;
        Ft = cosd( theta ) ;
        alpha = 1 ;
        
        F = FepF(Ft,alpha,Ib,Mxe,Vb,q) ;
        Isp = IspepF(Ft,alpha,mue,Mxe,Vb,q) ;
        fprintf( 'The thrust is %f N \n' , F ) ;
        fprintf( 'The specific impulse is %f seconds \n' , Isp )
        
        mdot = F/(Isp*9.81) ;
        tb = massp/mdot ;
        It = tb*F ;
        fprintf( 'The total operation time is %f hours \n' , tb/3600 ) 
        fprintf( 'The total impulse is %f kN*s \n' , It*1e-3 )
    end
end