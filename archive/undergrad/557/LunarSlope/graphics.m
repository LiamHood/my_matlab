clc ; close all ; clc ;
points{1} = [1,1,0] ;
points{2} = [1,0,0] ;
points{3} = [1,-1,0] ;
points{4} = [0,1,0] ;
points{5} = [0,0,0] ;
points{6} = [0,-1,0] ;
points{7} = [-1,1,0] ;
points{8} = [-1,0,0] ;
points{9} = [-1,-1,0] ;
linex = [ -1 , 0 , 1 ] ;
liney = [ -1 , 0 , 1 ] ;
lines{1} = [ linex(1) , liney(1) ; linex(3) , liney(3) ] ;
lines{2} = [ linex(3) , liney(1) ; linex(1) , liney(3) ] ;
lines{3} = [ linex(2) , liney(1) ; linex(2) , liney(3) ] ;
lines{4} = [ linex(1) , liney(2) ; linex(3) , liney(2) ] ;

figure
axis( [ -2 , 2 , -2 , 2 ] )
hold on
for ii = 1:9
    plot( points{ii}(1),points{ii}(2) , '*k')
end
% plot( 0 , 0 , '*b' )
for ii = 1:4
    plot( lines{ii}(:,1),lines{ii}(:,2) , 'r')
end
plot( 0 , 0 , 'ob' )