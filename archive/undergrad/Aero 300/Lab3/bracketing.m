function [bracket] = bracketing( f , n , g , h )
% i := 0
% a0 := initial guess
% while sign(f(ai)) =/= sign(f(a_i+1))
%       a_i+1 := a_i + h 
%       i := i + 1
% end
% a:a_i ; 
% b = a_i+1 ;


bracket = zeros(n,2) ; %preallocates

for jj = 1:n %repeats for assumed number of roots
    ii = 1 ; 
    a(1) = g(jj) ; %first boundary is the initial guess
    while sign( f( a(ii) ) ) == sign( f( a(ii) - h ) ) %continues as long as each bracket is on
        a(ii+1) = a(ii) + h ;
        ii = ii + 1 ;
    end
 
    bracket( jj , 1 ) = a( ii - 1 ) ;
    bracket( jj , 2 ) = a( ii ) ;
end

end
