function [w, D] = kohonen(w, x, alpha)
    n = size(w,2);
    D = zeros(1,n);
    for ii = 1:n
        D(ii) = sum((w(:,ii)-x).^2);
    end
    [~,ind] = min(D);
    w(:,ind) = w(:,ind) + alpha*(x-w(:,ind));
end