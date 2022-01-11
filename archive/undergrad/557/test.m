clear ; close all ; clc ;
RotMat = @(u,ang) [ cos(ang) + u(1)^2*(1-cos(ang)) , u(1)*u(2)*(1-cos(ang)) - u(3)*sin(ang) , u(1)*u(3)*(1-cos(ang)) + u(2)*sin(ang) ; ...
                    u(2)*u(1)*(1-cos(ang)) + u(3)*sin(ang) , cos(ang) + u(2)^2*(1-cos(ang)) , u(2)*u(3)*(1-cos(ang)) - u(1)*sin(ang) ; ...
                    u(3)*u(1)*(1-cos(ang)) - u(2)*sin(ang) , u(3)*u(2)*(1-cos(ang)) + u(1)*sin(ang) , cos(ang) + u(3)^2*(1-cos(ang)) ] ;

Rotx = @(ang) [1 0 0; 0 cos(ang) -sin(ang); 0 sin(ang) cos(ang)] ;
Roty = @(ang) [cos(ang) 0 sin(ang); 0 1 0; -sin(ang) 0 cos(ang)] ;
Rotz = @(ang) [cos(ang) -sin(ang) 0; sin(ang) cos(ang) 0; 0 0 1] ;
AxisRot = @(R,ang) [ R(3,2)-R(2,3) ; R(1,3)-R(3,1) ; R(2,1)-R(1,2) ]/(2*sin(ang)) ;
AngRot = @(R) acos( ( trace(R) - 1 )/2 ) ;
d2r = pi/180 ;
tspan = [ 0 , 1 ]*86400 ;
tol = 1e-10 ;

mu = 4905 ;
inci = 30*d2r ;
ecc = 1e-9 ;
raani = 10*d2r ;
omega = 0 ;
theta = 0 ;
a = 8000 ;
h = sqrt( mu*a*(1-ecc^2) ) ;
[ri,vi] = coes2state([h,inci,ecc,raani,omega,theta],mu) ;
[tiv, riv, viv] = TwoBody(tspan, ri, vi, mu, tol) ;

[req,veq] = coes2state([h,0,ecc,0,0,0],mu) ;
[teqv, reqv, veqv] = TwoBody(tspan, req, veq, mu, tol) ;

incf = 60*d2r ;
dinc = incf - inci ;
raanf = 50*d2r ;
draan = raanf - raani ;
[rf,vf] = coes2state([h,incf,ecc,raanf,omega,theta],mu) ;
[tfv, rfv, vfv] = TwoBody(tspan, rf, vf, mu, tol) ;



ba = acos( cos( inci )*cos( incf ) + sin( inci )*sin( incf )*cos( draan ) );
aoli = acos( ( sin( inci )*cos( draan ) - cos( ba )*sin( inci ) )/( sin( ba )*cos( inci ) ) )+pi/2 ;
[rint,vint] = coes2state([h,inci,ecc,raani,omega,aoli],mu) ;
urequ = [ cos(raani) ; sin(raani) ; 0 ] ;
% Rinc = RotMat(urequ,dinc) ;
% Rraan = Rotz(draan) ;
% Rmat = Rraan*Rinc ;
% urint = Rmat*urequ ;
% ang = AngRot(Rmat) ;
% uaxis = AxisRot(Rmat,ang) ;
% urequ = [ cos(raani) ; sin(raani) ; 0 ] ;
% guaxis = Rmat*uaxis ;
% rrint = Rmat*rint ;
% for ii = 1:length(tiv)
%     rrot(:,ii) = Rmat*riv(:,ii) ;
% end

inc3 = ba ;
[r3,v3] = coes2state([h,inc3,ecc,0,omega,theta],mu) ;
[t3v, r3v, v3v] = TwoBody(tspan, r3, v3, mu, tol) ;

RotFinc = Rotx(inci) ;
RotFraan = Rotz(raani) ;
hv = cross(riv(:,1),viv(:,1))/norm(cross(riv(:,1),viv(:,1))) ;
uinter = acos( ( sin( incf )*cos( draan ) - cos( ba )*sin( inci ) )/( sin( ba )*cos( inci ) ) ) ;
RotFu = RotMat(hv,uinter) ;
RmatF = RotFu*RotFraan*RotFinc ;

for ii = 1:length(t3v)
    r4v(:,ii) = RmatF*r3v(:,ii) ;
end
for ii = 1:length(teqv)
    reqr(:,ii) = RmatF*reqv(:,ii) ;
end

[sx,sy,sz] = sphere ;
figure
surf( sx*1e3 , sy*1e3 , sz*1e3 ) 
axis equal
hold on
plot3(reqv(1,:),reqv(2,:),reqv(3,:),'k')                % equitorial
plot3(riv(1,:),riv(2,:),riv(3,:),'b')                   % initial orbit
plot3(r3v(1,:),r3v(2,:),r3v(3,:),'m')                   % inc only
plot3(rfv(1,:),rfv(2,:),rfv(3,:),'r')                   % final orbit
% plot3(rrot(1,:),rrot(2,:),rrot(3,:),'g-')                % rotated initial orbit
plot3(r4v(1,:),r4v(2,:),r4v(3,:),'g')                   % rotation of 3
plot3(reqr(1,:),reqr(2,:),reqr(3,:),'c')                % rotation of equator
% plot3(rint(1,:),rint(2,:),rint(3,:),'*m')               % intersection point
% plot3(rrint(1,:),rrint(2,:),rrint(3,:),'*m')            % rotated intersection point
% plot3(a*urequ(1,:),a*urequ(2,:),a*urequ(3,:),'*b')      % ascending node
% plot3(a*urint(1,:),a*urint(2,:),a*urint(3,:),'*b')      % rotation of ascending node
% plot3(a*uaxis(1,:),a*uaxis(2,:),a*uaxis(3,:),'*g')      % axis of rotation
% plot3(a*guaxis(1,:),a*guaxis(2,:),a*guaxis(3,:),'*g')   % rotation of axis of rotation
