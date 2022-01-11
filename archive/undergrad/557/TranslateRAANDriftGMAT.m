clear ; close all ; clc ;
% h , inc , ecc , RAAN , omega , theta , a , rp , ra
name = 'YoD_FC_9000' ;
nameO = [ name , '.txt' ];
nameS = [ name , '.mat' ] ;
SatState = readmatrix( nameO ) ;
days = SatState(:,1) ;
len = length( days ) ;
hP = zeros( len , 1 ) ;
incP = zeros( len , 1 ) ;
eccP = zeros( len , 1 ) ;
raanP = zeros( len , 1 ) ;
aopP = zeros( len , 1 ) ;
thetaP = zeros( len , 1 ) ;
% raanP(1,1) = 360 ;
% raanP(1,5) = 360 ;
% raanP(1,9) = 360 ;
% raanP(1,13) = 360 ;
% raanP(1,17) = 360 ;
% raanP(1,21) = 360 ;
% raanP(1,25) = 360 ;
% raanP(1,29) = 360 ;

for jj = 0:7
    for ii = 1:4
        hP(:,ii+4*jj) = SatState(:,ii+1+25*jj) ;
        incP(:,ii+4*jj) = SatState(:,ii+5+25*jj) ;
        eccP(:,ii+4*jj) = SatState(:,ii+9+25*jj) ;
        for kk = 1:len
            if SatState(kk,ii+13+25*jj) == 0
                SatState(kk,ii+13+25*jj) = 360 ;
            end
            if ( SatState(1,ii+13+25*jj) - SatState(kk,ii+13+25*jj) ) >= 0
                raanP(kk,ii+4*jj) = SatState(kk,ii+13+25*jj) ;
            else
                raanP(kk,ii+4*jj) = SatState(kk,ii+13+25*jj) - 360 ;
            end
        end
        aopP(:,ii+4*jj) = SatState(:,ii+17+25*jj) ;
        thetaP(:,ii+4*jj) = SatState(:,ii+21+25*jj) ;
    end
end
% for ii = 1:4 %5-8
%     hP(:,ii+4) = SatState(:,ii+25) ;
%     incP(:,ii+4) = SatState(:,ii+29) ;
%     eccP(:,ii+4) = SatState(:,ii+33) ;
%     raanP(:,ii+4) = SatState(:,ii+37) ;
%     aopP(:,ii+4) = SatState(:,ii+41) ;
%     thetaP(:,ii+4) = SatState(:,ii+45) ;
% end
% for ii = 1:4 %5-8
%     hP(:,ii+8) = SatState(:,ii+25) ;
%     incP(:,ii+8) = SatState(:,ii+29) ;
%     eccP(:,ii+8) = SatState(:,ii+33) ;
%     raanP(:,ii+8) = SatState(:,ii+37) ;
%     aopP(:,ii+8) = SatState(:,ii+41) ;
%     thetaP(:,ii+8) = SatState(:,ii+45) ;
% end

save( nameS ) ;