%% Project
clear ; close all ; clc ;
lat = 35.3 ;
long = -120.66 ;
alt = 105.8 ;
mu = 398600 ;
r2d = 180/pi ;
coesName = [ "h" ; "inclination" ; "eccentricity" ; "RAAN" ; "argument of perigee" ; "true anomaly" ; "semi-major axis" ; "radius of periapsis" ; "radius of apoapsis" ] ;
        coesUnits = [ "km^2/s" ; "degrees" ; "  " ; "degrees" ; "degrees" ; "degrees" ; "km" ; "km" ; "km" ] ;
        
        
%% Part 1
JDbeg = juliandate( [ 2020 2 11 2 0 0 ] ) ;
JDend = juliandate( [ 2020 2 11 3 0 0 ] ) ;
tstep = 10 ;
[ r0 , v0 , epoch0 ] = TLE2state( 'ISSTLE_test.txt' ) ;
% [ r0 , v0 , epoch0 ] = TLE2state( 'Lacrosse5rocketTLE.txt' ) ;
COESinitial = state2coes_display( [ r0 ; v0 ] , mu ) ;
day = epoch0 - 31 ;
JDtle = juliandate( [ 2020 , 2 , day ] ) ;
[ t , r , v ] = Encke( 10 , [JDtle,JDbeg]*(3600*24) , r0 , v0 , mu , 'gravity,J2,J3' , 1 , 1 ) ;
% [ t , r , v ] = Encke( 10 , [0,JDbeg-JDtle]*(3600*24) , r0 , v0 , mu , '0' , 1 , 1 ) ;
[ JD , rho , az , el , VIS ] = PredictLookAnglesODE( r(:,end) , v(:,end) , lat , long , alt , JDbeg , JDend , tstep ) ;
% [ JD , rho , az , el , VIS ] = PredictLookAngles( r0 , v0 , lat , long , alt , JDtle , JDend , tstep ) ;
time = ( JD - JDbeg )*24*60 ;
figure
plot( time , rho )
% axis( [ 0 , 60 , 0 , 20000 ] )
figure
plot( time , az*(180/pi) )
% axis( [ 0 , 60 , 0 , 360 ] )
figure
plot( time , el*(180/pi) )
% axis( [ 0 , 60 , 0 , 90 ] )
figure 
plot( az*(180/pi) , el*(180/pi) )
% axis( [ 180 , 360 , 0 , 90 ] )

%% Part 3, 4, 5
ra = [ 30.381 ; 65.134 ; 99.976 ] ;
dec = [ 23.525 ; 0.774 ; -30.44 ] ;
UTC = [ 2013 , 03 , 25 , 03 , 10 , 30.032 ; 2013 , 03 , 25 , 03 , 15 , 20.612 ; 2013 , 03 , 25 , 03 , 20 , 32.777 ] ;
[ r2vec , v2vec , r1vec , r3vec ] = AnglesToRV( ra , dec , UTC , lat , long , alt ) ;
COES = state2coes_display( [ r2vec ; v2vec ] , mu ) ;
coesOut = [ coesName , COES' , coesUnits ] ;
for ii = 1:3
    JDsp(ii) = juliandate( UTC(ii,:) ) ;
end
tf = ( JDsp(2) - JDsp(1) )*86400 ;
tm = 1 ;
[ ~ , v2izzo ] = Lambert_Izzo( r1vec , r2vec , tf , tm , mu ) ;
[ ~ , v2bat ] = Lambert_Battin( r1vec , r2vec , tf , tm , mu ) ;
[ ~ , v2gau ] = Lambert_Gauss( r1vec , r2vec , tf , tm , mu ) ;
[ ~ , v2uv ] = Lambert_UV( r1vec , r2vec , tf , tm , mu ) ;
difi = norm( v2izzo - v2vec ) ;
difb = norm( v2bat - v2vec ) ;
difg = norm( v2gau - v2vec ) ;
difu = norm( v2uv - v2vec ) ;

disp( coesOut ) ;
disp( 'All of the inclinations are near each other so I can''t use that to ' )
disp( 'determine the best fitting object. RAAN is only a few degrees ' )
disp( 'off of object 2 so fits the best. RAAN is plane defining so should ' )
disp( 'be most accurate. Prop all TLEs to same time and compare. ' )

radpt = [62.0650000000000;190.293000000000;208.218000000000;218.144000000000;111.941000000000;127.243000000000;149.566000000000;173.968000000000] ;
decdpt = [53.3710000000000;55.8450000000000;33.9850000000000;14.0240000000000;29.5280000000000;24.4710000000000;10.9960000000000;-10.3980000000000] ;
UTCdpt = [2013,5,11,3,37,38.4680000000000;2013,5,11,3,44,41.3660000000000;2013,5,11,3,46,48.5990000000000;2013,5,11,3,49,6.78900000000000;2013,5,11,5,45,23.6950000000000;2013,5,11,5,47,30.8670000000000;2013,5,11,5,50,19.1180000000000;2013,5,11,5,53,31.1630000000000];
radp{1} = [ radpt(1) ; radpt(2) ; radpt(4) ] ;
radp{2} = [ radpt(1) ; radpt(3) ; radpt(4) ] ;
radp{3} = [ radpt(5) ; radpt(6) ; radpt(7) ] ;
radp{4} = [ radpt(5) ; radpt(7) ; radpt(8) ] ;
decdp{1} = [ decdpt(1) ; decdpt(2) ; decdpt(4) ] ;
decdp{2} = [ decdpt(1) ; decdpt(3) ; decdpt(4) ] ;
decdp{3} = [ decdpt(5) ; decdpt(6) ; decdpt(7) ] ;
decdp{4} = [ decdpt(5) ; decdpt(7) ; decdpt(8) ] ;
UTCdp{1} = [ UTCdpt(1,:) ; UTCdpt(2,:) ; UTCdpt(4,:) ] ;
UTCdp{2} = [ UTCdpt(1,:) ; UTCdpt(3,:) ; UTCdpt(4,:) ] ;
UTCdp{3} = [ UTCdpt(5,:) ; UTCdpt(6,:) ; UTCdpt(7,:) ] ;
UTCdp{4} = [ UTCdpt(5,:) ; UTCdpt(7,:) ; UTCdpt(8,:) ] ;
JDend = juliandate( UTCdpt(8,:) ) ;
for ii = 1:4
    [ rdp{ii} , vdp{ii} , r1dp{ii} , r3dp{ii} ] = AnglesToRV( radp{ii} , decdp{ii} , UTCdp{ii} , lat , long , alt ) ;
    JD(ii) = juliandate( UTCdp{ii}(2,:) ) ;
    [ t , rfv{ii} , vfv{ii} ] = CowellJ2J3( [JD(ii),JDend+10]*(3600*24) , rdp{ii} , vdp{ii} , mu ) ;
    
    COESdp(ii,:) = state2coes_display( [ rfv{ii}(:,end) ; vfv{ii}(:,end) ] , mu ) ;
end
for ii = 1:9
    COESavg(ii) = mean( COESdp(:,ii) ) ;
end
    coesOut = [ coesName , COESavg' , coesUnits ] ;
    disp( coesOut ) ;

[ robj{1} , vobj{1} , epoch{1} ] = TLE2state( 'P1_dp1.txt' ) ;
[ robj{2} , vobj{2} , epoch{2} ] = TLE2state( 'P1_dp2.txt' ) ;
[ robj{3} , vobj{3} , epoch{3} ] = TLE2state( 'P1_dp3.txt' ) ;
for ii = 1:1
    epoch{ii} = juliandate( [ 2013 , 1 , 1 , 0 , 0 , 0 ] ) + epoch{ii} - 1 ;
    disp( JDend - epoch{ii} )
    [ t , rofv , vofv ] = CowellJ2J3( [epoch{ii},JDend+10]*(3600*24) , robj{ii} , vobj{ii} , mu ) ;
%     [ t , rofv , vofv ] = TwoBody( [epoch{ii},JDend]*(3600*24) , robj{ii} , vobj{ii} , mu , 1e-8 ) ;
    disp([ 'Object ' , num2str( ii ) ] );
    COESdpo{ii} = state2coes_display( [ rofv(:,end) ; vofv(:,end) ] , mu ) ;
    coesOuto = [ coesName , COESdpo{ii}' , coesUnits ] ;
    disp(coesOuto) ;
end


