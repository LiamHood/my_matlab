clear; close all; clc;

problem1()
problem1a()

problem2()

%% Functions
function problem1()
    fprintf("\nProblem 1:\n")
    A = [-1.10188, 0.90528, -.00212;
            4.0639, -.77013, -.169190;
            0, 0, 10];
    B = [0, 0; 0, 1; 10, 0];
    [vec, poles] = eig(A);
   fprintf("Open Loop Poles are \\lamda = \n\t%f + %fi \n\t%f + %fi \n\t%f + %fi\n", ...
        [real(poles(1,1)), imag(poles(1,1)), real(poles(2,2)), imag(poles(2,2)), ...
        real(poles(3,3)), imag(poles(3,3))])

    cont_matrix = [B, A*B, A^2*B];
    fprintf("Controllable since controllability matrix is rank %i and A is rank %i\n", [rank(cont_matrix), rank(A)])
    
    lam1 = -2+2i;
    des_lambda = [ lam1, conj(lam1), -15];
    I3 = eye(3);
    D = diag([1, 1, 1]);
    vec1 = [0.20 + .35i, -.98 + .07i, 0]';
    des_vec = [vec1, conj(vec1), [0; 0; 1]];
    uv = @(lambda,vec) pinv([lambda.*I3 - A, B; D, zeros(3,2)])*[0;0;0;vec];
    for ii = 1:3
        uvii = uv(des_lambda(ii), des_vec(:,ii));
        vd(:,ii) = uvii(1:3);
        ud(:,ii) = uvii(4:5);
    end
    K = ud*inv(vd);
    fprintf("[\t%f + %fi,\t%f + %fi,\t%f + %fi\n",real(K(1,:)), imag(K(1,:)))
    fprintf(" \t%f + %fi,\t%f + %fi,\t%f + %fi]\n",real(K(2,:)), imag(K(2,:)))
    [vec, poles] = eig(A-B*K);
    fprintf("Closed Loop Poles are \\lamda = \n\t%f + %fi \n\t%f + %fi \n\t%f + %fi\n", ...
        [real(poles(1,1)), imag(poles(1,1)), real(poles(2,2)), imag(poles(2,2)), ...
        real(poles(3,3)), imag(poles(3,3))])
end

function problem1a()
    fprintf("\nProblem 1:\n")
    A = [-1.10188, 0.90528, -.00212;
            4.0639, -.77013, -.169190;
            0, 0, 10];
    B = [0, 0; 0, 1; 10, 0];
    sys = ss(A, B, eye(3), 0);
    [vec, poles] = eig(A);
   fprintf("Open Loop Poles are \\lamda = \n\t%f + %fi \n\t%f + %fi \n\t%f + %fi\n", ...
        [real(poles(1,1)), imag(poles(1,1)), real(poles(2,2)), imag(poles(2,2)), ...
        real(poles(3,3)), imag(poles(3,3))])

    cont_matrix = [B, A*B, A^2*B];
    fprintf("Controllable since controllability matrix is rank %i and A is rank %i\n", [rank(cont_matrix), rank(A)])
    
    lam1 = -2+2i;
    des_lambda = [ lam1, conj(lam1), -15];
    I3 = eye(3);
    D = diag([1, 1, 1]);
    vec1 = [0.20 + .35i, -.98 + .07i, 0]';
    des_vec = [vec1, conj(vec1), [0; 0; 1]];
    uv = @(lambda,vec) pinv([lambda.*I3 - A, B; D, zeros(3,2)])*[0;0;0;vec];
    for ii = 1:3
        uvii = uv(des_lambda(ii), des_vec(:,ii));
        vd(:,ii) = uvii(1:3);
        ud(:,ii) = uvii(4:5);
    end
    K = ud*inv(vd);
    fprintf("[\t%f + %fi,\t%f + %fi,\t%f + %fi\n",real(K(1,:)), imag(K(1,:)))
    fprintf(" \t%f + %fi,\t%f + %fi,\t%f + %fi]\n",real(K(2,:)), imag(K(2,:)))
    [vec, poles] = eig(A-B*K);
    fprintf("Closed Loop Poles are \\lamda = \n\t%f + %fi \n\t%f + %fi \n\t%f + %fi\n", ...
        [real(poles(1,1)), imag(poles(1,1)), real(poles(2,2)), imag(poles(2,2)), ...
        real(poles(3,3)), imag(poles(3,3))])
end

function problem2()
    fprintf("\nProblem 2:\n")
    A = [-1.10188, 0.90528, -.00212;
           4.0639, -.77013, -.169190;
           0,      0,       10];
    B = [0,  0;
         0,  1; 
         10, 0];
    C = [1, 0, 0;
         0, 1, 0];
    [vec, poles] = eig(A);
    fprintf("Open Loop Poles are \\lamda = \n\t%f + %fi \n\t%f + %fi \n\t%f + %fi\n", ...
        [real(poles(1,1)), imag(poles(1,1)), real(poles(2,2)), imag(poles(2,2)), ...
        real(poles(3,3)), imag(poles(3,3))])

    cont_matrix = [B, A*B, A^2*B];
    fprintf("Controllable since controllability matrix is rank %i and A is rank %i\n", [rank(cont_matrix), rank(A)])
    
    lam1 = -2+2i;
    des_lambda = [ lam1, conj(lam1)];
    I3 = eye(3);
    D = [1, 0, 0;
         0, 1, 0];
    O22 = zeros(2,2);
    vec1 = [.2+.35i, -.98+.07i]';
    des_vec = [vec1, conj(vec1)];
    uv = @(lambda,vec) pinv([lambda.*I3 - A, B; D, O22])*[0;0;0;vec];
    for ii = 1:2
        uvii = uv(des_lambda(ii), des_vec(:,ii));
        vd(:,ii) = uvii(1:3);
        ud(:,ii) = uvii(4:5);
    end
    K = ud*inv(C*vd);
    fprintf("[\t%f + %ii,\t%f + %ii\n",real(K(1,:)), imag(K(1,:)))
    fprintf(" \t%f + %ii,\t%f + %ii]\n",real(K(2,:)), imag(K(2,:)))
    [vec, poles] = eig(A-B*K*C);
    fprintf("Closed Loop Poles are \\lamda = \n\t%f + %fi \n\t%f + %fi \n\t%f + %fi\n", ...
        [real(poles(1,1)), imag(poles(1,1)), real(poles(2,2)), imag(poles(2,2)), ...
        real(poles(3,3)), imag(poles(3,3))])
end
