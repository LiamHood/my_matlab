clear ; close all ; clc ;

% silicon
    % Stopping powers
        % electrons
        spes = 1.531 ;
        % protons
        spps = 1.754e2 ;
        % alpha
        spas = 1.296e3 ;
    rhoskg = 2300 ; % density of silicon in kg/m^3 
    rhos = rhoskg*(1e3)/(1e2)^3 ; % g/cm^3
    % Penetration Depth
        Res = 1/(rhos*spes) ;
        Rps = 1/(rhos*spps) ;
        Ras = 1/(rhos*spas) ;
disp( 'a' )
disp( 'The stopping power of alpha particles is largest so it will do the ')
disp( 'most damage but electrons will reach the furthest ')
disp( 'Penetration Depth in Silicon' )
disp([ 'Electrons  ' , num2str(Res) , ' meters' ])
disp([ 'Protons  ' , num2str(Rps) , ' meters' ])
disp([ 'Alpha  ' , num2str(Ras) , ' meters' ])
disp( ' ' )
% Aluminum
    % Stopping powers
        % electrons
        spea = 1.486 ;
        % protons
        sppa = 1.720e2 ;
        % alpha
        spaa = 1.226e3 ;
    rhoakg = 2700 ; % density of aluminum in kg/m^3 
    rhoa = rhoakg*(1e3)/(1e2)^3 ; % g/cm^3
    % Penetration Depth
        Rea = 1/(rhoa*spea) ;
        Rpa = 1/(rhoa*sppa) ;
        Raa = 1/(rhoa*spaa) ;
disp( 'b' )
disp( 'The stopping power of alpha particles is largest so it will do the ')
disp( 'most damage but electrons will reach the furthest ')
disp( 'Penetration Depth in Aluminum' )
disp([ 'Electrons  ' , num2str(Rea) , ' meters' ])
disp([ 'Protons  ' , num2str(Rpa) , ' meters' ])
disp([ 'Alpha  ' , num2str(Raa) , ' meters' ])
disp( ' ' )
% Graphite
    % Stopping powers
        % electrons
        speg = 1.627 ;
        % protons
        sppg = 2.297e2 ;
        % alpha
        spag = 1.893e3 ;
    rhogkg = 2150 ; % density of silicon in kg/m^3 
    rhog = rhogkg*(1e3)/(1e2)^3 ; % g/cm^3
    % Penetration Depth
        Reg = 1/(rhog*speg) ;
        Rpg = 1/(rhog*sppg) ;
        Rag = 1/(rhog*spag) ;
disp( 'b' )
disp( 'The stopping power of alpha particles is largest so it will do the ')
disp( 'most damage but electrons will reach the furthest ')
disp( 'Penetration Depth in Graphite' )
disp([ 'Electrons  ' , num2str(Reg) , ' meters' ])
disp([ 'Protons  ' , num2str(Rpg) , ' meters' ])
disp([ 'Alpha  ' , num2str(Rag) , ' meters' ])

disp(' ' )
disp( 'The aluminum functions best as a radiation shield as the radiation' )
disp( 'penetrates the least far in it' )


