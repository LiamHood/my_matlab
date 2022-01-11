clear ; close all ; clc ;
% Load Data
load( 'CH_data.mat' )

% Generate Satellite Info
P = 4 ;
a = 7500 ;
i = 55 ;
F = 2 ;
S = T/P ;
PU = 360/T ;
for kk = 1:T
	Plane = kk - floor( kk/P )*P ; 
	SatInPlane = floor( kk/P ) ;
	Constellation(kk).A = a ;
	Constellation(kk).E = 1e-9 ;
	Constellation(kk).I = i ;
	Constellation(kk).W = 90 ;
					
	raanInitial(kk) = S*PU*kk ;
	Constellation(kk).RAAN = mod( raanInitial(kk) , 360 ) ;
	Constellation(kk).TA = PU*F*Plane + PU*P*SatInPlane ;
end

% Change Large Matrix to 3-D Matrix
satGround = zeros( 288 , T , n ) ;
for kk = 1:n
    for ii = 1:T 
        for jj = 1:288
            offset = (kk-1)*288 ;
            satGround(jj,ii,kk) = satSeen( jj + offset , ii ) ;
        end
    end
end

% Find time of first and second pass in a day
for ii = 1:n            % check all points
    [ I , J ] =  find( satGround(:,:,ii) == 1 ) ; 
    for jj = 1:T        % check all satellites
        start = min( find( J == jj ) ) ; % start of first pass
        kk = 0 ;
        visTime(ii,jj) = 0 ;
        if start < length( J )
            dif = abs( I(kk+start+1) - I(kk+start) ) ;
        else
            dif = 2 ;
        end
        while dif < 2
            visTime(ii,jj) = 300 + visTime(ii,jj) ;
            if kk + start + 1 <= length(I)
                dif = abs( I(kk+start+1) - I(kk+start) ) ;
            else 
                dif = 2 ;
            end
            kk = kk + 1 ;
        end
        if kk ~= 0 
            if kk + start + 1 < length( I ) %check the second pass to see if it is longer
                visTime2 = 0 ;
                if J(kk+start) == jj
                    start2 = start + kk ;
                    kk = 0 ;
                    dif = abs( I(kk+start2+1) - I(kk+start2) ) ;
                    while dif < 2
                        visTime2 = 300 + visTime2 ;
                        if kk + start2 + 1 <= length(I)
                            dif = abs( I(kk+start2+1) - I(kk+start2) ) ;
                        else
                            dif = 2 ;
                        end
                        kk = kk + 1 ;
                    end
                end
            else
                visTime2 = 0 ;
            end
        else
            visTime2 = 0 ;
        end
        if visTime2 > visTime(ii,jj)
            visTime(ii,jj) = visTime2 ; % replace first pass if second pass is longer
        end
    end
end

for ii = 1:n
    COV_avg(ii) = mean( COV(:,ii) ) ;
    COV_min(ii) = min( COV(:,ii) ) ;
end
for ii = 1:n
    timeSatPass(ii) = prctile( visTime(ii,:) , 75 ) ;
end

% ii = points, jj = satellites
points = 1:n ;
sats = 1:T ;


figure
plot( timeSatPass/3600 , '.' )
xlabel( 'Point' )
ylabel( 'Pass Time [hours]' )

figure
plot( COV_avg.*2 , '.' )
xlabel( 'Point' )
ylabel( 'Average Satellites in View' )

for ii = 1:length( COV_min )
    if COV_min(ii) == 0
        COV_min(ii) = .5 ;
    end
end
figure
plot( COV_min.*2 , '.' )
xlabel( 'Point' )
ylabel( 'Minimum Satellites in View' )

lowPassTime = prctile( timeSatPass , 8 ) ;
lowCOV = prctile( COV_avg , 8 ) ;

lowCoveragePoints = find( timeSatPass <= lowPassTime & COV_avg <= lowCOV ) ;

fprintf( 'Points with low pass times and low coverage are: \n' )
for ii = 1:length( lowCoveragePoints )
    fprintf( '%i \n' , lowCoveragePoints(ii) )
end


while (1)
    fprintf('Please input point of Interest\n');

    point = input('? ');

    if isempty( point )
    else
        if point == 'q'
            break ;
        else
            
            strTitle = [ 'Point: ' , num2str( point ) , ' Lat: ' , num2str( pointLoc( point , 1 ) ) , ' Long: ' , num2str( pointLoc( point , 2 ) ) ] ;
            med = median( visTime(point,:)/3600 ) ;
            sev = prctile( visTime(point,:)/3600,75 ) ;
            twe = prctile( visTime(point,:)/3600,25 ) ;
            medLine = [med,med] ;
            sevLine = [sev,sev] ;
            tweLine = [twe,twe] ;

            figure
            hold on
            plot( [0,T] , sevLine , 'b')
            plot( [0,T] , medLine , 'r')
            plot( [0,T] , tweLine , 'g')
            plot( sats(:) , visTime(point,:)./3600 , '.k' ) 
            legend( '75th percentile' , '50th percentile' , '25th percentile' , 'Pass Time' , 'Location' , 'East')
            hold off
            xlabel( 'Satellite #' )
            ylabel( 'Pass Time [hours]' )
            title( strTitle ) ;
        end
    end
end

