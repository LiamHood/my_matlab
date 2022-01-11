function ap = PerturbedAccelerations( forces , r , v , A , m , JDo , t )    
            if contains( forces , "drag" )
                if contains( forces , "NRLMSISE" )
                    [ epochnew , yds ] = epochUpdate( epoch , t ) ;
                    apd = DragNRLMSISE( r , v , A , m , epochnew , yds ) ;
                else
                    apd = DragAcceleration( r , v , A , m ) ;
                end
            else
                apd = zeros(3,1) ;
            end
            if contains( forces , "gravity" )
                if contains( forces , "J2" )
                    aph = J2accel( r ) ; 
                end
                if contains( forces , "J3" )
                    aph = aph + J3accel( r ) ; 
                end
            else 
                aph = zeros(3,1) ; 
            end
            if contains( forces , "nbody" )
                apn = ThreeBody( t , r , JDo , 'Sun' ) ;
                
            else
                apn = zeros(3,1) ;
            end
            if contains( forces , "srp" )
                aps = SRPacceleration( r , v , A , m , JDo , t , 1 ) ;
            else
                aps = zeros(3,1) ;
            end
            ap = apd + aph + apn + aps ;

end