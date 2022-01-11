load( 'Cylinder_Lab_25_Set4.mat' )
P_25_4 = P ;
load( 'Cylinder_Lab_25_Set3.mat' )
P_25_3 = P ;
Re_25_3 = Re ;
load( 'Cylinder_Lab_25_Set2.mat' )
P_25_2 = P ;
Re_25_2 = Re ;
load( 'Cylinder_Lab_25_Set1.mat' )
P_25_1 = P ;
Re_25_1 = Re ;
P_25_avg_full = ( P_25_1+P_25_2+P_25_3 ) / 3 ;
P_25_avg_full(:,2) = (P_25_avg_full(:,1)+P_25_avg_full(:,3))./2 ;
for ii = 3:26
    P_25_avg(ii-2) = mean(P_25_avg_full(:,ii));
end
theta = 15 ;
for ii = 2:24
    theta(ii) = theta(ii-1) + 15 ;
end
P_25_avg(2) = (P_25_avg(1)+P_25_avg(3))/2 ;
plot( theta , P_25_avg )
theta
