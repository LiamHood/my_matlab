function [f] = fun(p,x)
%turns p vector of coefecients into function
f = 0;
    lop = length(p) ;
    for kk = 1:lop
    f = f + x^(lop-kk) * p( kk ) ;
    end

end

