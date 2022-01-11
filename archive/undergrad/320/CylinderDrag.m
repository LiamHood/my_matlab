load( 'Cy_500_1.mat' )
d2r = 180/pi ;
theta = 15*d2r ;
for ii = 1:24
    f(ii) = - p(ii) * r * sin( theta ) * theta ;
end