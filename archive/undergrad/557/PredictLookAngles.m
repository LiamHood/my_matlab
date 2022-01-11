function [ JD , rho , az , el , VIS ] = PredictLookAngles( r0 , v0 , lat , long , alt , JDbeg , JDend , tstep ) 
    mu = 398600 ;
    r = r0 ;
    v = v0 ;
    JD = JDbeg ;
    ii = 2 ;
    while JD( ii-1 ) < JDend
%         [ r(:,ii) , v(:,ii) ] = NewStateUV( r(:,ii-1) , v(:,ii-1) , tstep , mu );
        [ t , r , v ] = TwoBody( tspan , r0 , v0 , mu , tol ) ;
        time = julian2date( JD(ii-1) ) ;
        time(6) = time(6) + tstep ;
        JD(ii) = juliandate( time ) ;
        ii = ii + 1 ;
    end
    for ii = 1:length( JD )   
        UTC = julian2date( JD(ii) ) ;
        [ rho(ii,1) , az(ii,1) , el(ii,1) , ~ , ~ , ~ , VIS(ii,1) ] = RAZEL( r(:,ii) , v(:,ii) , UTC , lat , long , alt ) ;
    end
end