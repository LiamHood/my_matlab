clear ; close all ; clc ;
% visibility = readmatrix( 'SeeEarth.txt' ) ;
% visibility = [ visibility(:,1) , visibility(:,3:53) ] ;
% save( 'Visibility.mat' )
load( 'Visibility.mat' )
time = visibility(:,1) ;
ii = 1 ;
while time(ii) <= 1
    for jj = 1:51 
        if visibility(ii,jj) == 0 
            satsee(ii,jj) = .0035*24*60 ;
        else
            satsee(ii,jj) = 0 ;
        end
    end
    ii = ii + 1 ;
end
for ii = 1:51 
    covered(1,ii) = sum( satsee(:,ii) ) ;
end
ii = find( time == 1 ) ;
while time(ii) <= 2
    for jj = 1:51 
        if visibility(ii,jj) == 0 
            satsee(ii,jj) = .0035*24*60 ;
        else
            satsee(ii,jj) = 0 ;
        end
    end
    ii = ii + 1 ;
end
for ii = 1:51 
    covered(2,ii) = sum( satsee(:,ii) ) ;
end
ii = find( time == 2 ) ;
while time(ii) <= 3
    for jj = 1:51 
        if visibility(ii,jj) == 0 
            satsee(ii,jj) = .0035*24*60 ;
        else
            satsee(ii,jj) = 0 ;
        end
    end
    ii = ii + 1 ;
end
for ii = 1:51 
    covered(3,ii) = sum( satsee(:,ii) ) ;
end
ii = find( time == 3 ) ;
while time(ii) <= 4
    for jj = 1:51 
        if visibility(ii,jj) == 0 
            satsee(ii,jj) = .0035*24*60 ;
        else
            satsee(ii,jj) = 0 ;
        end
    end
    ii = ii + 1 ;
end
for ii = 1:51 
    covered(4,ii) = sum( satsee(:,ii) ) ;
end
ii = find( time == 4 ) ;
while time(ii) <= 5
    for jj = 1:51 
        if visibility(ii,jj) == 0 
            satsee(ii,jj) = .0035*24*60 ;
        else
            satsee(ii,jj) = 0 ;
        end
    end
    ii = ii + 1 ;
end
for ii = 1:51 
    covered(5,ii) = sum( satsee(:,ii) ) ;
end