function dy = euler_fun(t, y, J, T_c)

dy = zeros(10,1);

w = y(1:3);
E = y(4:6);
q_vec = y(7:9);
q4 = y(10);

dy(1:3,1) = J\(T_c - cross(w, J*w));

dy(4:6,1) = LBI(E(1), E(2), E(3),313)*w;

dy(7:9) = 1/2*(q4*w - cross(w, q_vec));
dy(10) = -1/2*(w'*q_vec);