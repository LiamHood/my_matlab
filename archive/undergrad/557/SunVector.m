function rs = SunVector( JD )
    Tutc = ( JD - 2451545.0 )/364525 ;
    lambdaMS = 280.460 + 36000.771*Tutc ;
    Ttbd = Tutc ;
    MS = 357.5291092 + 35999.05034*Tutc ;
    lambdaEcl = lambdaMS + 1.914666471*sin( MS ) + 0.019994643*sin( 2*MS ) ;
    rsmag = ( 1.000140612 - 0.016708617*cos( MS ) - 0.000139589*cos( 2*MS ) )*149597870.7 ;
    eps = 23.439291 - 0.0130042*Ttbd ;
    rs = [ rsmag*cos( lambdaEcl ) ; rsmag*cos( eps )*sin( lambdaEcl ) ; rsmag*sin( eps )*sin( lambdaEcl ) ] ;
end