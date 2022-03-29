function [y, a] = maxnet(x, epsilon)
    n = length(x);
    a = zeros(1,n);
    y = zeros(1,n);
    for ii = 1:n
        a(ii) = x(ii) - epsilon*(sum(x)-x(ii));
        if a(ii) > 0
            y(ii) = a(ii);
        end
    end
end