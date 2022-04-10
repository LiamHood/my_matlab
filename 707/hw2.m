clear; close all; clc;

problem1()

problem2()

problem3()

problem4()


%% Functions

function problem1()
    fprintf("\nProblem 1:\n")
    A = [-1.10188, 0.90528, -.00212;
           4.0639, -.77013, -.169190;
            0,           0,      -10];
    B = [0, 0;
         0, 1;
         10, 0];
    C = eye(3);

    lam1 = -2+2i;
    lam2 = conj(lam1);
    lam3 = -15;
    des_poles = [lam1, lam2, lam3];

    vec1 = [0.20 - .35i, -.98 - .07i, 0]';
    vec2 = conj(vec1);
    vec3 = [0; 0; 1];
    des_vec = [vec1, vec2, vec3];
    
    EigStructAssign(A, B, C, des_poles, des_vec);
    
end

function problem2()
    fprintf("\nProblem 2:\n")
    A = [-1.10188, 0.90528, -.00212;
           4.0639, -.77013, -.169190;
           0,      0,       -10];
    B = [0,  0;
         0,  1; 
         10, 0];
    C = [1, 0, 0;
         0, 1, 0];
    
    lam1 = -2+2i;
    lam2 = conj(lam1);
    des_poles = [lam1, lam2];

    vec1 = [0.20 - .35i, -.98 - .07i]';
    vec2 = conj(vec1);
    des_vec = [vec1, vec2];
    
    EigStructAssign(A, B, C, des_poles, des_vec);
    fprintf("The acuator pole is moved from -15 to 9.5 now that we are not \n")
    fprintf("able to directly assign the actuator pole\n\n")
end

function problem3()
    fprintf("\nProblem 3:\n")
    A = [-1.10188, 0.90528, -.00212;
           4.0639, -.77013, -.169190;
           0,      0,       -10];
    B = [0,  0, 10]';
    C = [1, 0, 0;
         0, 1, 0];
    
    lam1 = -2+2i;
    lam2 = conj(lam1);
    des_poles = [lam1, lam2];

    vec1 = [0.20 - .35i, -.98 - .07i]';
    vec2 = conj(vec1);
    des_vec = [vec1, vec2];
    
    EigStructAssign(A, B, C, des_poles, des_vec);
    fprintf("The actuator pole has moved increasingly positive to 12.15\n\n")
end

function problem4()
fprintf("\nProblem 4:\n")
    A = [0, 1;
        -2, -3];
    poles = eig(A);
    fprintf("Open Loop Poles are \\lamda = \n")
    disp(poles)
    fprintf("The poles are negative so the open loop system is stable\n")
    B = [0; 2];
    Ackermann(A, B, [-3, -5]);
end

function K = Ackermann(A, B, lambda)
    n = length(A);

    syms s
    funs = 1;
    for lam = lambda
        funs = funs*(s-lam);
    end
    fdes_s = collect(funs);
    coef = sym2poly(fdes_s);
    fdes_A = coef(1).*A*A + coef(2).*A + coef(3).*eye(n);

    Con = [];
    for ii = 1:n
        Con = [Con, A^(ii-1)*B];
    end

    leftv = [zeros(1,n-1), 1];
    K = leftv*inv(Con)*fdes_A;
    fprintf("Control Gains are K = \n")
    disp(K)
end

function K = EigStructAssign(A, B, C, des_lambda, des_vec)
    poles = eig(A);
    fprintf("Open Loop Poles are \\lamda = \n")
    disp(poles)
    
    cont_matrix = [B, A*B, A^2*B];
    fprintf("Controllable since controllability matrix is rank %i and A is rank %i\n\n", [rank(cont_matrix), rank(A)])

    n = length(A);
    m = size(B,2);
    p = length(des_lambda);

    I = eye(n);
    D = diag(ones(1,p));
    O = zeros(p,(n+m)-p);
    top = zeros(n,1);
    uv = @(lambda,vec) pinv([lambda.*I - A, B; D, O])*[top;vec];

    for ii = 1:p
        uvii = uv(des_lambda(ii), des_vec(:,ii));
        vd(:,ii) = uvii(1:n);
        ud(:,ii) = uvii(n+1:end);
    end
    K = ud*inv(C*vd);
    fprintf("Control Gains are K = \n")
    disp(K)
    poles = eig(A-B*K*C);
    fprintf("Closed Loop Poles are \\lamda = \n")
    disp(poles)
    fprintf("Desired Eigen Vectors are v = \n")
    disp(des_vec)
    fprintf("Achievable Eigen Vectors are v = \n")
    disp(vd)
end