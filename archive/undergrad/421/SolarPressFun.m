function Torque = SolarPressFun( Cbg , r , s , Areas , Positions , Norms )

r = r/norm(r) ;

sundirection = dot( r , s ) ;

theta = acosd(sundirection) ;

lam = asind( 6378 / norm(r) ) ; % https://ocw.tudelft.nl/wp-content/uploads/AE2104-Orbital-Mechanics-Slides_10.pdf

eclipseang = 180-lam/2 ;

if theta <= eclipseang 

    %At 1 AU, solar pressure is 

    ps = 2*( 1361/299792458 ) ; %Pascals

    %The normal vectors on the body frame are pulled in

    %Use rotation matrix to rotate from ECI to Bodyframe

    sb = Cbg*s ; %rotates the sun vector into the body frame

    %Find the amount of area on each side, exposed to the sun, use dot products

    Forces = 0 ; %preallocates force summation

    Top = [0;0;0] ; %Preallocates Cps Parts 

    Bot = 0 ; 

    for pp = 1:length(Areas)

        dotted = (Norms(:,pp)'*s) ; %Dot product of norm and sun

        if dotted > 0 

            Forces = Forces + -ps*sb*dotted*Areas(pp) ;

        end 

        % Need Center of Pressure Now


        if dotted > 0 

            Top = Top + Positions(:,pp) * dotted * Areas(pp) ; % Sums the position vector times the dotted value and area

            Bot = Bot + dotted * Areas(pp) ; %Sums the n dot s * area

        end

    end

    Cps = Top/Bot ; %Defines the center of pressure

    Torque = cross(Cps,Forces) ; %Cross product of force and distance

else

    Torque = [ 0 ; 0 ; 0 ] ;

end

end