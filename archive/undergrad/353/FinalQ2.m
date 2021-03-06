%% Final Question 2
% Liam Hood
clear; close all; clc;

% Given
rad = 1 ;
z = 400 ;
Tplas = 3500 ;
denplas = 1e12 ;
sacol = pi ; % area ions are collected over
vsc = sqrt( 398600 / ( z + 6378 ) )*1000 ;

kb = 1.3806e-23 ;
qe = 1.60217662e-19 ;
me = 9.10938356e-31 ;

disp( 'Assumptions ' )
disp('     Ignore Ise, Isi, Ibse ')
disp('     eclipse Iph=0 ')
disp('     uniform s/c so Ic = 0 ')
disp('     ions cannot impact wake side')
disp('     quasi-neutral so ne = ni ')
disp('     single ions only so charge of electron is same a charge of ion ')
disp( ' ' )
V = linspace( -5 , -1 , 1e2 ) ;
Ie = .25*qe*denplas*(8*kb*Tplas/(pi*me))*sacol*4*exp(qe*V/(kb*Tplas)) ;
Ii = qe*denplas*vsc*sacol ;
Ii = ones( 100 , 1 )*Ii ;

figure
subplot( 2 , 1 , 1 )
semilogy( V , Ie )
title( 'Electron Current' )
xlabel( 'Potential (V)' )
ylabel( 'Current' )
subplot( 2 , 1 , 2 )
plot( V , Ii )
title( 'Ion Current' )
xlabel( 'Potential (V)' )
ylabel( 'Current' )

Vfloat = (kb*Tplas)/(qe) * log(vsc*sqrt((pi*me)/(8*kb*Tplas))) ;
disp([ 'The voltage will drift to ' , num2str( Vfloat ) , ' Volts' ])
