function [ JD , rho , az , el , VIS ] = PredictLookAnglesODE( r0 , v0 , lat , long , alt , JDbeg , JDend , tstep ) 
    mu = 398600 ;
    r = r0 ;
    v = v0 ;
    JD = JDbeg ;
    ii = 2 ;
    tspan = [ 0 , JDend-JDbeg ]*(24*60*60) ;
    tol = 1e-8 ;
    [ t , r , v ] = TwoBody( tspan , r0 , v0 , mu , tol ) ;
    JD = JDbeg + t/86400 ;
    for ii = 1:length( JD )   
        UTC = julian2date( JD(ii) ) ;
        [ rho(ii,1) , az(ii,1) , el(ii,1) , ~ , ~ , ~ , VIS(ii,1) ] = RAZEL( r(:,ii) , v(:,ii) , UTC , lat , long , alt ) ;
    end
end