function [ RAANdot , omegadot ] = oblatenessPerturbation( mu , J2 , R , ecc , a , inc )
% Inputs and outputs in degrees
RAANdotr = -((3/2)*((sqrt(mu)*J2*R^2)/((1-ecc^2)^2*a^(7/2))))*cosd(inc) ;
omegadotr = RAANdotr*(((5/2)*sind(inc)^2-2)/cosd(inc)) ;
RAANdot = RAANdotr*(180/pi) ;
omegadot = omegadotr*(180/pi) ;
end