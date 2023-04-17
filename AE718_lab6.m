clear; close all; clc;

% y1 = 586;
% y3 = -119;
% y2 = 120;
% y4 = 1425;
% 
% syms I t1 t2 t3 t4
% 
% % eq1 = t1 + t2 + t3 +t4 == 360;
% eq2 = y1 == I*cosd(t1);
% eq3 = y2 == I*cosd(t2);
% eq4 = y3 == I*cosd(t3);
% eq5 = y4 == I*cosd(t4);
% 
% answer = solve([eq2,eq3,eq4,eq5],[I,t1,t2,t3,t4])
% at(1) = answer.t1;
% at(2) = answer.t2;
% at(3) = answer.t3;
% at(4) = answer.t4;
% for ii = 1:4
%     while at(ii) > 360
%         at(ii) = at(ii) - 360;
%     end
%     while at(ii) < -360
%         at(ii) = at(ii) + 360;
%     end
% end
% at
syms a b c intensity_correction

light = [a;b;c];
yaw1 = [1;0;0];
yaw2 = [0;1;0];
yaw3 = [-1;0;0];
yaw4 = [0;-1;0];
top = [0;0;-1];

ry1 = [586;296;175;175;84];
ry3 = [-119;-146;-122;-122;-101];
ry2 = [120;160;140;140;110];
ry4 = [1425;960;712;712;348];
rt = [216;89;24;24;8];

for ii = 1:5
    equ1(ii) = a*intensity_correction == ry1(ii);
    equ2(ii) = b*intensity_correction == ry2(ii);
    equ3(ii) = -a*intensity_correction == ry3(ii);
    equ4(ii) = -b*intensity_correction == ry4(ii);
    equ5(ii) = -c*intensity_correction == rt(ii);
    equ6(ii) = norm(light) == 1;
    solution(ii) = vpasolve([equ1(ii),equ4(ii),equ5(ii),equ6(ii)],[a,b,c,intensity_correction]);
    light_direction(:,ii) = double([solution(ii).a;solution(ii).b;solution(ii).c]);
end
light_direction

% syms a b c intensity_correction
% 
% light = [a;b;c];
% yaw1 = [1;0;0];
% yaw2 = [0;1;0];
% yaw3 = [-1;0;0];
% yaw4 = [0;-1;0];
% top = [0;0;-1];
% 
% ry1 = 22.5;
% ry4 = 48.04;
% rt = 9.62;
% % for ii = 1
% %     equ1(ii) = dot(light,yaw1)*intensity_correction == ry1(ii);
% %     equ2(ii) = dot(light,yaw2)*intensity_correction == ry2(ii);
% %     equ3(ii) = dot(light,yaw3)*intensity_correction == ry3(ii);
% %     equ4(ii) = dot(light,yaw4)*intensity_correction == ry4(ii);
% %     equ5(ii) = dot(light,top)*intensity_correction == rt(ii);
% %     solution(ii) = solve([equ1(ii),equ2(ii),equ3(ii),equ4(ii),equ5(ii)],[a,b,c]);
% % end
% for ii = 1
%     equ6 = norm(light) == 1;
%     equ1(ii) = a*intensity_correction == ry1(ii);
%     equ4(ii) = -b*intensity_correction == ry4(ii);
%     equ5(ii) = -c*intensity_correction == rt(ii);
%     solution(ii) = vpasolve([equ1(ii),equ4(ii),equ5(ii),equ6(ii)],[a,b,c,intensity_correction]);
% end
% light_direction = [solution.a;solution.b;solution.c]
% 
