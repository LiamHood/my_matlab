function rootB = bisecting( a , b , TOL , f )
% given [ a b ] such that f(a)*f(b) < 0 will find the roots of 
% while ( b - a ) / 2 > TOL
    % c = ( b + a ) / 2
    % if f(c) is 0 
        %stop
    % end
    % if sign of f(c)*f(b) is positive ;
        % b = c
    % else
        % a = c
    % end
% end
for ii = 1:length(a)
    while ( b(ii) - a(ii) ) / 2 > TOL %Continue running as long as half the difference between b and a is greater than the tolerance 
        c(ii) = ( b(ii) + a(ii) ) / 2 ; % c is halway between a and b
        if f(c(ii)) == 0 %if c is the root the function ends
            stop
        end
        if f(c(ii))*f(b(ii)) > 0 %if the f(c) and f(b) are on the same side of the x axis
            b(ii) = c(ii) ; % new b is now c
        else %otherwise
            a(ii) = c(ii) ; % the lower bound is now c
        end
    end
    rootB(ii) = c(ii) ;
end
end


