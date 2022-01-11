function f = fund( pf , x , d )

lp = length(pf) ;
fm = zeros( lp ) ; %allocates space for the function matrix
    for kk = 1:d 
         for mm = 1:lp
            fm(kk,mm) = ( x^( lp - mm ) * pf( kk , mm ) ); %Evaluates each cell by multiplying the value by the appropriate value of x
         end
    end
    
    f = sum(fm') ; %adds all the terms of the function matrix

end