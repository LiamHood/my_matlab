function [ v , dd ] = orbitsHW1( v , r  )
dx = v(1) ;
dy = v(2) ;
dz = v(3) ;
ddx = -(mu_e/norm(r)^3)*r(1) ;
ddy = -(mu_e/norm(r)^3)*r(2) ;
ddz = -(mu_e/norm(r)^3)*r(3) ;
v = [ dx dy dz ] ;
dd = [ ddx ddy ddz ] ;

end