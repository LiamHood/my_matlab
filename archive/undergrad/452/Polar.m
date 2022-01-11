%--------------------------------------------------------------------------
% Az  : azimuth of vector
% Elev: altitude of vector
% r   : norm of vector
%--------------------------------------------------------------------------
function Vec = Polar(Az, Elev, r)

if (nargin == 2)
    r = 1;
end

Vec = zeros(3,1);

Vec(1,1) = r * cos(Az) * cos(Elev);
Vec(2,1) = r * sin(Az) * cos(Elev);
Vec(3,1) = r * sin(Elev);

