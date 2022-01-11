clear
tic
p = [ 1 , 1.0142 , -19.3629 , 15.8398 ];


function rootL = allroot( p , TOL , start , n )
%Finds all roots greater than starting value
%Step with increasing size until finding root then bisect to find root
    %each time through the step size increases 
    %

ii = 1 ; 
jj = 1 ;
atr = start ;
ss = size( atr );
while ss(2) < n %repeats if less roots than assumed
    h = .01 ; %initial step size
    
    %adds a column to fill in every time it repeats
    if ii > 1
        if ii <= n
            atr = [ atr , zeros(ss(1),1) ] ;
        end
    end
    
    %starts the next iteration of root where the last root was found
    if ii > 1 
        atr(1,ii) = atr( jj , ii - 1 ) ;
    end
    
    jj = 1 ;
        %Bracketing but with increasing step size
        while sign( fun( p , atr(jj,ii) ) ) == sign( fun( p , atr(jj,ii) - h ) ) %Stops when values are found to be on either side of a root
            ss = size(atr) ;
            
            if  jj + 1 == ss(1)
                atr = [ atr ; zeros( 1 , ss(2) ) ] ; %adds a row if there is not one
            end
            
            atr(jj + 1,ii) = atr(jj,ii) + h ; %Increase boundary by step size
            jj = jj + 1 ;
            h = h + .1 ;
            
        end
        
        if fun( p , atr( jj , ii ) ) == 0 %stop if root is already found
            
        else
            b = atr( jj , ii ) ; 
            a = atr( jj - 1 , ii ) ; 
            
            while ( b - a ) / 2 > TOL %Continue running as long as half the difference between b and a is greater than the tolerance 
                c = ( b + a ) / 2 ; % c is halway between a and b
                
                if fun(p,c) == 0 %if c is the root the function ends
                    stop
                end
                
                if fun(p,c)*fun(p,b) > 0 %if the f(c) and f(b) are on the same side of the x axis
                    b = c ; % new b is now c
                else %otherwise
                    a = c ; % the lower bound is now c
                end
                
            end

            root(ii) = c ;
        
            
        end
        

ii = ii + 1 ;
ss = size( atr );
end

    %alltheroots = x(:,ii) 


toc

