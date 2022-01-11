clear ; close all ; clc ;
d2r = pi/180 ;
inci = 30*d2r ;
incf = 60*d2r ;
raani = 30*d2r ;
raanf = 90*d2r ;
draan = raanf - raani ;
dinc = incf - inci ;

RotMat = @(u,ang) [ cos(ang) + u(1)^2*(1-cos(ang)) , u(1)*u(2)*(1-cos(ang)) - u(3)*sin(ang) , u(1)*u(3)*(1-cos(ang)) + u(2)*sin(ang) ; ...
                    u(2)*u(1)*(1-cos(ang)) + u(3)*sin(ang) , cos(ang) + u(2)^2*(1-cos(ang)) , u(2)*u(3)*(1-cos(ang)) - u(1)*sin(ang) ; ...
                    u(3)*u(1)*(1-cos(ang)) - u(2)*sin(ang) , u(3)*u(2)*(1-cos(ang)) + u(1)*sin(ang) , cos(ang) + u(3)^2*(1-cos(ang)) ] ;

% Rotx = @(ang) [1 0 0; 0 cos(ang) -sin(ang); 0 sin(ang) cos(ang)] ;
% Roty = @(ang) [cos(ang) 0 sin(ang); 0 1 0; -sin(ang) 0 cos(ang)] ;
Rotz = @(ang) [cos(ang) -sin(ang) 0; sin(ang) cos(ang) 0; 0 0 1] ;
AxisRot = @(R) [ R(3,2)-R(2,3) ; R(1,3)-R(3,1) ; R(2,1)-R(2,1) ] ;
AngRot = @(R) acos( ( trace(R) - 1 )/2 ) ;

ui = [ cos(raani) ; sin(raani) ; 0 ] ;
Rinc = RotMat(ui,dinc) ;
Rraan = Rotz(draan) ;
Rmat = Rinc*Rraan ;
% Rmat = Rraan*Rinc ;
totdp = AngRot(Rmat)/d2r ;


% ba = acos( cos( inci )*cos( incf ) + sin( inci )*sin( incf )*cos( draan ) );
% ui = acos( ( sin( incf )*cos( draan ) - cos( ba )*sin( inci ) )/( sin( ba )*cos( inci ) ) ) ;
% uf = acos( ( cos( inci )*sin( incf ) - sin( inci )*cos( incf )*cos( draan ) )/sin( ba ) ) ;
% uid = ui/d2r ;
% ufd = uf/d2r ;
% du = uf - ui ;
% tli = (ui + raani)/d2r ;
% tlf = (uf + raanf)/d2r ;
% ali = ( inci*sin( ui ) )/d2r ;
% alf = ( incf*sin( uf ) )/d2r ;
