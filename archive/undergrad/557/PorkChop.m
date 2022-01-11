clear ; close all ; clc ;
mu = 1.327124e11 ;
JDbegE = 2453528.0 ;
JDendE = 2453682.0 ;
JDbegM = 2453706.0 ;
JDendM = 2454156.0 ;
JDE = JDbegE:1:JDendE ;
JDM = JDbegM:1:JDendM ;
n = length( JDE ) ;
m = length( JDM ) ;
for ii = 1:n
    for jj = 1:m 
%         datiE = julian2date( JDE(ii) ) ;
%         datiM = julian2date( JDM(jj) ) ;
        [~ , rEv(:,ii) , vEv(:,ii) , ~ ] = AERO557planetcoe_and_sv( 3 , JDE(ii) ) ;
        [~ , rMv(:,jj) , vMv(:,jj) , ~ ] = AERO557planetcoe_and_sv( 4 , JDM(jj) ) ;
        rE = rEv(:,ii) ;
        vE = vEv(:,ii) ;
        rM = rMv(:,jj) ;
        vM = vMv(:,jj) ;
        dt(ii,jj) = ( JDM(jj) - JDE(ii) )*86400 ;
        [ vEt , vMt ] = Lambert_Izzo( rE , rM , dt(ii,jj) , 1 , mu ) ;
        dvE = norm( vEt - vE ) ;
        dvM = norm( vMt - vM ) ;
        dv(ii,jj) = dvE + dvM ;
        vinf(ii,jj) = dvM ;
        c3(ii,jj) = dvE^2 ;
        [ vEtL , vMtL ] = Lambert_Izzo( rE , rM , dt(ii,jj) , -1 , mu ) ;
        dvEL = norm( vEtL - vE ) ;
        dvML = norm( vMtL - vM ) ;
        dvL(ii,jj) = dvEL + dvML ;
        vinfL(ii,jj) = dvML ;
        c3L(ii,jj) = dvEL^2 ;
        if dvL(ii,jj) < dv(ii,jj)
            dv(ii,jj) = dvL(ii,jj) ;
            vinf(ii,jj) = vinfL(ii,jj) ;
            c3(ii,jj) = c3L(ii,jj) ;
        end
            
    end
end
levdv = [ 1:.5:9 , 12:3:21 ] ;
levvinf = [ 4.5 , 6:3:21 ] ;
levc3 = [ 0:1:18 , 20 , 30:20:100 ] ;
levdt = 0:50:600 ;

figure
contour( JDE - JDbegE , JDM - JDbegM , dv' , levdv , 'r' ,  'ShowText' , 'on'  )
hold on
contour( JDE - JDbegE , JDM - JDbegM , vinf' , levvinf , 'b' , 'ShowText' , 'on' )
contour( JDE - JDbegE , JDM - JDbegM , c3' , levc3 , 'g' , 'ShowText' , 'on' )
contour( JDE - JDbegE , JDM - JDbegM , dt'./86400 , levdt , 'k' , 'ShowText' , 'on' )

figure
contour( JDE - JDbegE , JDM - JDbegM , dv' , levdv , 'r' ,  'ShowText' , 'on'  )
figure
contour( JDE - JDbegE , JDM - JDbegM , vinf' , levvinf , 'b' , 'ShowText' , 'on' )
figure
contour( JDE - JDbegE , JDM - JDbegM , c3' , levc3 , 'g' , 'ShowText' , 'on' )
figure
contour( JDE - JDbegE , JDM - JDbegM , dt'./86400 , levdt , 'k' , 'ShowText' , 'on' )
% figure
% hold on
% plot3( rEv(1,:) , rEv(2,:) , rEv(3,:) )
% plot3( rMv(1,:) , rMv(2,:) , rMv(3,:) )