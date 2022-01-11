%      % FIX THIS
     for ii = 1:length(t)
%      if t(ii) >= 60
%          tmin = floor( t(ii)/60 ) ;
%          t(ii) = t(ii) - tmin(ii)*60 ;
%          if tmin(ii) >= 60
%              thour(ii) = floor( tmin(ii)/60 ) ;
%              tmin(ii) = tmin(ii) - thour(ii)*60 ;
%          else
%              thour(ii) = 0 ;
%          end
%      else
%          tmin(ii) = 0 ;
%          thour(ii) = 0 ;
%      end
%      if ( thour(ii) + UTC(4) ) >= 24 
%          thour(ii) = thour(ii) - 24 ;
%          UTC(3) = UTC(3) + 1 ;
%      end
%      UTCnow = UTC + [ 0 0 0 thour(ii) tmin(ii) t(ii) ] ;

    JDay = Day0 + t(ii)/(24*60*60) ; % Find Julian day of current moment
    JYear = JDay/365.2422 ; % Find Julian year 
    [ rerels(:,ii) , ~ ] = planetary_state(3,JDay) ; %Find position vector of earth relative to the sun
     rads = norm(rerels) ;
     end
     figure
     plot3(rerels(1,:),rerels(2,:),rerels(3,:))
     figure 
     plot( rads )
     norm(rerels)