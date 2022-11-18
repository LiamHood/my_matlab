syms x t p(x) u(x) b a

eq1 = int(u(x)*diff(p(x), x)*diff(u(x), x), x, a, b)