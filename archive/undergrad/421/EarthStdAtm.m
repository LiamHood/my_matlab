function [rho] = EarthStdAtm(alt)

% Chart Values Provided by the website below

% http://toughsf.blogspot.com/2017/09/low-earth-orbit-atmospheric-scoops.html

Altitudes = [50 60 70 80 100 200 300 400 500 600 700 800 900 1000] ; % Chart altitudes in km

Densities = [1.03E-3 3.10E-4 8.28E-5 1.85E-5 5.36E-7 3.13E-10 2.40E-11 3.38E-12 6.21E-13 1.39E-13 4.03E-14 1.66E-14 9.11E-15 5.85E-15] ; %Densities for given altitudes in kg/m^3

% Use linear interpolation

Pos = find(Altitudes>alt) ;

if alt>1000

    warning('AHHHHHHH, altitude too high!')

end

Dens1 = Densities(Pos(1)) ; %Finds the density that is right below input

Dens2 = Densities(Pos(1)+1) ; %Finds density above input altitude

AltDiff = Altitudes(Pos(1)+1)-Altitudes(Pos(1)) ; % Finds Difference in altitudes

rho = Dens1 + (Dens2-Dens1)/AltDiff ; %Linearly interpolates the rho at a given altitude


end
