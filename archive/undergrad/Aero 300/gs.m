%Check for convergence
%Reiterate as long as error is less than preceeding error
    %Run for each row
        %Calculate first sum
        %calculate second sum
        %calculate new x guess
    %Check to ensure continuing convergence
    %reiterate
    
clear

a = [ 8 1 -2 -3 ; -1 -10 2 5 ; 1 -6 12 -3 ; -3 2 3 -9 ] ;
b = [ 1 ; 2 ; 3 ; 4] ;
TOL = .01 ;

n = length(a);
    
%Set starting values
    x = ones( n , 1 ) ;
    
    kk = 1 ;
    convergence = 1;
    err(kk) = [ 10^3err ] ;
    ax = 0 ;
    %Convergence test 
    for mm = 1:n %runs following for every row
        rowsum = 0 ;
        for oo = 1:n
            rowsum = rowsum + abs(a(mm,oo)); %Finds the sum of each row
        end
            if mm == oo %Runs for all elements on main diagonal
                if abs(a(mm,oo)) <  rowsum - abs(a(mm,oo)) %displays warning if the diagonal value is less that the rowsum
                    warning( 'May not converge' )
                end
            end
        
    end
    
    while convergence > 0
        for ii = 1:n
            sum1 = 0 ;
            sum2 = 0 ;
            for jj = 1 : (ii-1) %First sum for each row
                sum1 = sum1 + a(ii,jj)*x(jj,kk+1) ;
            end
            for jj = (ii+1) : n %Second sum for each row
                sum2 = sum2 + a(ii,jj)*x(jj,kk) ;
            end
            x(ii,kk+1) = ( 1/a(ii,ii) ) * ( b(ii) - sum1 - sum2 ) ;
        end
        
        err(kk+1) = norm( a*x(:,kk) ) ;
        convergence = err(kk) - err(kk+1) ;
        
%         if kk > 1
%             err(kk+1) = norm( x( 1:n , kk ) - x( 1:n , kk - 1 ) ) ;
%         else
%             err(kk+1) = err(kk);
%         end
        x( 1:n , kk+2 ) = zeros( n , 1 ) ;
        kk = kk + 1 ;
    end
    