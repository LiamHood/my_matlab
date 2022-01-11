clear; close all; clc;
% 
%    data = readmatrix( 'Lab2Data.xlsx' ) ;
%    save( 'Lab2Data.mat' ) ;
%     
load( 'Lab2Data.mat' ) ;
load( 'Langmuir_practice.txt' )
trial{1} = data( : , 1:3 ) ;
trial{2} = data( : , 5:7 ) ;
trial{3} = data( : , 9:12 ) ;

    figure
    hold on 
for ii = 1:3 
    I(:,ii) = trial{ii}(:,3) ;
    V(:,ii) = trial{ii}(:,1) ;
    index0 = find( V(:,ii) == 0 ) ;
    for jj = 1:index0 
        I(jj,ii) = -I(jj,ii) ;
    end
    [ ioffset , index ] = min(I(:,ii)) ;
    I(:,ii) = I(:,ii) - ioffset ;
    Vfloat{ii} = V( index,ii ) ;

%     e = 1.60217662e-19 ;
%     me = 9.1938356e-31;
%     mi = 6.65e-26 ;
%     k = 1.38064852e-23 ;
%     Te{ii} = e*(V(index0+1,ii)-V(index0-1,ii))/(k*log(I(index0+1,ii)/I(index0-1,ii))) ; % kelvin
%     Ti{ii} = Te{ii} ;
%     Aprobe = 1 ;
%     veth{ii} = sqrt( 8*k*Te{ii}/me )*Aprobe;
%     vith{ii} = sqrt( k*Te{ii}/mi )*Aprobe;
%     n{ii} = ioffset./(.25*e*vith{ii}*Aprobe) ;
%     j{ii} = n{ii}*e*veth{ii} ;

    slope = (I(index0+1,ii)/I(index0-1,ii))/(V(index0+1,ii)-V(index0-1,ii)) ;
%     midline = slope.*V(:,ii) - I(index0,ii) ;


    plot( V(:,ii) , I(:,ii) , '.')
end
xlabel( 'Potential (mV)' )
ylabel( 'Current (A)' )

for ii = 1:length(I) 
    Iavg(ii) = mean( I(ii,:) ) ;
    Vavg(ii) = mean( V(ii,:) ) ;
end

    Vfloat = Vavg( index ) ;

    e = 1.60217662e-19 ;
    me = 9.1938356e-31;
    mi = 6.65e-26 ;
    k = 1.38064852e-23 ;
    eps0 = 8.85418782e-12 ;
    Te = e*(V(index0+1)-V(index0-1))/(k*log(I(index0+1)/I(index0-1))) ; % kelvin
    Ti = Te ;
    Aprobe = 1 ;
    veth = sqrt( 8*k*Te/me )*Aprobe;
    vith = sqrt( k*Te/mi )*Aprobe;
    n = -ioffset./(.25*e*vith*Aprobe) ;
    j = n*e*veth ;
    debye = sqrt( eps0*k*Te/(n*(e^2)) ) ;
    pp = (4/3)*pi*debye^3*n ;
    slope = (Iavg(index0+1)/Iavg(index0-1))/(Vavg(index0+1)-Vavg(index0-1)) ;
    Ies = (1/4)*e*n*veth*Aprobe ;    
    PlasPot = ( Ies + slope*Vavg(index0) - Iavg(index0) )/slope ;
    freq = 8.98*sqrt(n) ;
    
 plot( Vavg , Iavg , '.b' )
 legend( 'Trial 1' , 'Trial 2' , 'Trial 3' , 'Average' , 'Location' , 'northwest' )

 hold off
 
figure
hold on
 plot( Vavg , Iavg , '.b')
xlabel( 'Potential (mV)' )
ylabel( 'Current (A)' )


load( 'Lab2Data_G3.mat' )
V3 = data(:,1) ;
I3 = data(:,2)*(1e-3)/220 ;
    [ ioffset3 , index ] = min(I3) ;
    I3 = I3 - ioffset3 ;
plot( V3 , I3 , '.' )

    Vfloat3 = V3( index ) ;

    Te3 = e*(V3(index0+1)-V3(index0-1))/(k*log(I3(index0+1)/I3(index0-1))) ; % kelvin
    Ti3 = Te3 ;
    Aprobe = 1 ;
    veth3 = sqrt( 8*k*Te3/me )*Aprobe;
    vith3 = sqrt( k*Te3/mi )*Aprobe;
    n3 = -ioffset3./(.25*e*vith3*Aprobe) ;
    j3 = n*e*veth3 ;
    debye3 = sqrt( eps0*k*Te3/(n3*(e^2)) ) ;
    pp3 = (4/3)*pi*debye3^3*n3 ;
    slope3 = (I3(index0+1)/I3(index0-1))/(V3(index0+1)-V3(index0-1)) ;
    Ies3 = (1/4)*e*n3*veth3*Aprobe ;
    PlasPot3 = ( Ies3 + slope3*V3(index0) - I3(index0) )/slope3 ;
    freq3 = 8.98*sqrt(n3) ;
    
load( 'Lab2Data_G1.mat' )
V1 = data(:,1) ;
I1 = data(:,3) ;
    [ ioffset1 , index ] = min(I1) ;
    I1 = I1 - ioffset1 ;
plot( V1 , I1 , '.' )

    Vfloat1 = V1( index ) ;
    index0 = 32 ;
    Te1 = e*(V1(index0+1)-V1(index0-1))/(k*log(I1(index0+1)/I1(index0-1))) ; % kelvin
    Ti1 = Te1 ;
    veth1 = sqrt( 8*k*Te1/me )*Aprobe;
    vith1 = sqrt( k*Te1/mi )*Aprobe;
    n1 = -ioffset1./(.25*e*vith1*Aprobe) ;
    j1 = n*e*veth1 ;
    debye1 = sqrt( eps0*k*Te1/(n1*(e^2)) ) ;
    pp1 = (4/3)*pi*debye1^3*n1 ;
    slope1 = (I1(index0+1)/I1(index0-1))/(V1(index0+1)-V1(index0-1)) ;
    Ies1 = (1/4)*e*n1*veth1*Aprobe ;    
    PlasPot1 = ( Ies1 + slope1*V1(index0) - I1(index0) )/slope1 ;
    freq1 = 8.98*sqrt(n1) ;
    
legend( 'Our Data' , 'Group 3' , 'Group 1' , 'Location' , 'northwest' )
hold off
