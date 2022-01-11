clear ; clc ; close all ;
spaces = 1e2 ;
z = 450 ;
inc = 45 ;
S = 110 ;
k = 1 ;
d = logspace( -6 , 3 , spaces ) ;
H = ( 10.^exp( -( (log(d)-.78) ./ .637 ) .^2 ) ).^.5 ;
af1 = 10^( z/200 - S/140 -1.5 ) ;
af = af1 / ( af1 + 1 ) ;
incd = .990 ;
F1 = 1.22e-5*d.^-2.5 ;
F2 = 8.1e10*(d+700).^-6 ;
q = .02 ;
qprime = .04 ;
g1 = ( 1 + q )^23 * ( 1 + qprime )^( 2019 - 2011 ) ;
p = .05 ;
g2 = 1 + p*( 2019 - 1988 ) ;
Fod = H.*af.*incd.*( F1.*g1 + F2.*g2 ) ;

figure
loglog( d(1:spaces) , Fod )
xlabel( 'Diameter (cm)' )
ylabel( 'Flux (Particles/(m^2*yr))' )
title( 'Orbital Debris Flux at 450km at 45 degree inclination' )


