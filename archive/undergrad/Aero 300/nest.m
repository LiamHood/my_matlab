% Program 0.1 Nested Multiplication
% Evaluates the polynomial from nested form using Horner's Method
% From Numerial Analysis, Timothy Sauer

function y = nest(d, c, x, b)

 if nargin < 4
     b = zeros(d,1);
 end

y = c(d + 1);

for i = d:-1:1
    y = y.*(x-b(i))+c(i);
    %y = y.*x+c(i);
end
