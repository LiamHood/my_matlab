function [y, a] = mexican_hat(x, R1, w1, R2, w2)
    n = length(x);
    a = zeros(1,n);
    y = zeros(1,n);
    for ii = 1:n
        r1index = ii-R1:ii+R1;
        r2index = [ii-R2:ii-R1, ii+R1:ii+R2];
        a(ii) = w1*sum(x(r1index(r1index>0 & r1index<n))) ...
              + w2*sum(x(r2index(r2index>0 & r2index<n)));
        if a(ii) > 0
            y(ii) = a(ii);
        end
    end
end