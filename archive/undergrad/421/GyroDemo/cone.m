function h = cone(h, rx, ry, n, origin, direction, varargin)
%
% usage:
%
% [h] = cone(rx, ry, l)
%
% h = handle of cone graphics object
% rx = radius of cone in x direction (semi-axis)
% ry = radius of cone in y direction (semi-axis)
% l = a 3-D vector defining the axis of the cone.  The cone vertex is
% located at origin and opens in the direction of l.
%
% [h] = cone(..., 'property', 'value')
%
% ... = any other patch property pairs
%
% Eric A. Mehiel
% Cal Poly, SLO
% Aerospace Engineering Department

r = (0:.1:1);

[x, y, z] = cylinder(r, n);

x = rx*x;
y = ry*y;
%l = norm(direction);
l = h;
z = l*z;

if isempty(varargin)
    h = patch(surf2patch(x,y,z), 'edgeColor','none', 'faceColor', [.5 .5 .5]);
else
    h = patch(surf2patch(x,y,z), 'edgeColor','none', varargin{:});
end

translate(h, origin(1), origin(2), origin(3))

angle = acosd(dot([0 0 1], direction/norm(direction)));

if (angle == 180)
    axis = [1 0 0];
    rotate(h, axis, angle, origin);
elseif (angle ~= 0)
    axis = cross([0 0 1], direction);
    rotate(h, axis, angle, origin);
end

