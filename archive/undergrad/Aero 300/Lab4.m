clear
a = [ 8 1 -2 -3 ; -1 -10 2 5 ; 1 -6 12 -3 ; -3 2 3 -9 ] ;
b = [ 1 ; 2 ; 3 ] ;

x_j = jacobi( a , b )

%% Functions

%If not diagonally dominant throw a warning about possible non convergence
%before continuing

%Take an initial guess of x(1) to the solution x
    %for k greater than 1 continue
        %for ii from 1 to n 
            %x(ii) begins at 0
            %for jj from 1 to n
                %ax(ii) = ax(ii) plus a(ii,jj)*x(jj)
            %end
        %x(ii) = (b(ii) - ax(ii) + a(ii)*x(ii)
        %end
    %x^(k) = x 
    %check convergence
    %end
function [x_j] = jacobi( a , b )
    n = length(a);
    rowsum = 0 ;
    xi = zeros( n , 1 ) ;
    xf = zeros( n , 1 ) ;
    kk = 2 ;
    dif = [ 0 0 0 ] ;
    ax = 0 ;
    %Convergence test 
    for mm = 1:n
        for oo = 1:n
            rowsum = rowsum + a(mm,oo); %Finds the sum of each row
            if mm == oo
                if abs(a(mm,oo)) <  rowsum %displays warning if the diagonal value is less that the rowsum
                    warning( 'May not converge' )
                end
            end
        end
    end
    
    while dif(kk) < dif(kk-1)
        kk = kk + 1 ;
        xi = xf ;
        for ii = 1:n
            for jj = 1:n
                ax = ax + a(ii,jj)*xi(jj) ;
            end
            xf(ii) = ( 1/a(ii,ii) ) * ( b(ii) - ( ax - a(ii,ii)*xi(ii) ) ) ;
        end
        for mm = 1 : ii
            dif(kk-1) = abs( xi(ii) - xf(ii) ) ;
        end
        kk = kk + 1 ;
    end
    x_j = xf ;
end


function [x_gs] = gauss_seidel( a , b )

end
