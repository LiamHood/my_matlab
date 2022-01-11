function [ r2vec , v2vec , r1vec , r3vec ] = AnglesToRV( ra , dec , UTC , lat , long , alt )
        d2s = (3600*24) ; 
        qhat = zeros(3) ;
        rsite = zeros(3) ;
            for ii = 1:3 
                qhat(1,ii) = cosd( dec(ii) )*cosd( ra(ii) ) ;
                qhat(2,ii) = cosd( dec(ii) )*sind( ra(ii) ) ;
                qhat(3,ii) = sind( dec(ii) ) ;
                rsite(:,ii) = lla2eci( [ lat , long , alt ] , UTC(ii,:) )*1e-3 ;
            end
            jt = juliandate( UTC ) ;
            tau1 = ( jt(1) - jt(2) )*d2s ;
            tau3 = ( jt(3) - jt(2) )*d2s ;
        [ r2vec , v2vec , r1vec , r3vec ] = DoubleR( qhat(:,1) , qhat(:,2) , qhat(:,3) , rsite(:,1) , rsite(:,2) , rsite(:,3) , tau1 , tau3 ) ;
%         [ r2vec , v2vec , r1vec , r3vec ] = GaussExtended( qhat(:,1) , qhat(:,2) , qhat(:,3) , rsite(:,1) , rsite(:,2) , rsite(:,3) , tau1 , tau3 ) ;

end