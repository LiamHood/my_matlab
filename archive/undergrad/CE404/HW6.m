clear ; close all ; clc ;
stress = [-2625.40000000000;-2744.40000000000;-2863.40000000000;-2979.40000000000;-3095.30000000000;-3219.70000000000;-3344.10000000000;-3487;-3630;-3801.30000000000;-3972.60000000000;-4184;-4395.30000000000;-4661.10000000000;-4926.80000000000;-5266.40000000000;-5606;-6046.70000000000;-6487.40000000000;-7067.20000000000;-7647.10000000000;-8424.50000000000;-9201.80000000000;-10262;-11323;-12794;-14266;-16351;-18436;-21326;-24216;-28394;-32571;-36160;-39748;-41925;-44102;-44304] ;
zout = [0;0.263160000000000;0.526320000000000;0.789470000000000;1.05260000000000;1.31580000000000;1.57890000000000;1.84210000000000;2.10530000000000;2.36840000000000;2.63160000000000;2.89470000000000;3.15790000000000;3.42110000000000;3.68420000000000;3.94740000000000;4.21050000000000;4.47370000000000;4.73680000000000;5;5.26320000000000;5.52630000000000;5.78950000000000;6.05260000000000;6.31580000000000;6.57890000000000;6.84210000000000;7.10530000000000;7.36840000000000;7.63160000000000;7.89470000000000;8.15790000000000;8.42110000000000;8.68420000000000;8.94740000000000;9.21050000000000;9.47370000000000;9.73680000000000];
b = 1 ;
zdistnorm = -( 10 - zout )./b ;
stressnorm = stress./(-1e5) ;

plot( stressnorm , zdistnorm )
xlabel( 'Stress over applied Area Load' )
ylabel( 'Z distance over 1 edge of applied area' )

tmB = readmatrix( 'HW6mB.txt' ) ;
tmT = readmatrix( 'HW6mT.txt' ) ;
tmL = readmatrix( 'HW6mL.txt' ) ;
tmR = readmatrix( 'HW6mR.txt' ) ;
tmMv = readmatrix( 'HW6mMv.txt' ) ;
tmMh = readmatrix( 'HW6mMh.txt' ) ;

dmh(:,1) = tmB(:,1) ;
mmh(:,1) = tmB(:,2) ;
dmh(:,3) = tmT(:,1) ;
mmh(:,3) = tmT(:,2) ;
dmv(:,1) = tmL(:,1) ;
mmv(:,1) = tmL(:,2) ;
dmv(:,3) = tmR(:,1) ;
mmv(:,3) = tmR(:,2) ;
dmv(:,2) = tmMv(:,1) ;
mmv(:,2) = tmMv(:,2) ;
dmh(:,2) = tmMh(:,1) ;
mmh(:,2) = tmMh(:,2) ;

for ii = 1:3
    mmv(:,ii) = -mmv(:,ii)/(1*max( mmv(:,ii) )) ;
    mmh(:,ii) = -mmh(:,ii)/(1*max( mmh(:,ii) )) ;
end
% mmv = -mmv./max( max( mmh ) ) ;
% mmh = -mmh./max( max( mmh ) ) ;
mmh(:,2) = mmh(:,2) + 4.5 ;
mmh(:,3) = mmh(:,3) + 9 ;
mmv(:,2) = mmv(:,2) + 3 ;
mmv(:,3) = mmv(:,3) + 6 ;

figure
axis square
hold on
for ii = 1:3
    plot( mmh(:,ii) , dmh(:,ii) )
    plot( dmv(:,ii) , mmv(:,ii) )
end
title( 'Moment Distribution Abaqus' )


% 
% figure
% plot( dmB , mmB )
% xlabel( 'Distance Along Edge [m]' )
% ylabel( 'Moment [Nm/m]' )
% title( 'Bottom Edge' )
% 
% figure
% plot( mmL , dmL )
% ylabel( 'Distance Along Edge [m]' )
% xlabel( 'Moment [Nm/m]' )
% title( 'Left Edge' )
% 
% figure
% plot( dmT , mmT )
% xlabel( 'Distance Along Edge [m]' )
% ylabel( 'Moment [Nm/m]' )
% title( 'Top Edge' )
% 
% figure
% plot( mmR , dmR )
% ylabel( 'Distance Along Edge [m]' )
% xlabel( 'Moment [Nm/m]' )
% title( 'Right Edge' )
% 
% figure
% plot( dmMh , mmMh )
% xlabel( 'Distance Along Edge [m]' )
% ylabel( 'Moment [Nm/m]' )
% title( 'Middle Horizontal' )
% 
% figure
% plot( mmMv , dmMv )
% ylabel( 'Distance Along Edge [m]' )
% xlabel( 'Moment [Nm/m]' )
% title( 'Middle Vertical' )
% 
