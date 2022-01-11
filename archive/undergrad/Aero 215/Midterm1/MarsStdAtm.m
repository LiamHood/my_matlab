function [D0,P0,T0] = MarsStdAtm(h)
%input h in [m]
    T0 = 228.50; % [K]
    P0 = 610; % [Pa]
    D0 = 1.2250; % [kg/m^3]

    DT(1) = -1.9;
    DT(2) = 2.2;
    DT(3) = -1.4;
    DT(4) = -.5;
    DT(5) = -3.85;
    DT = DT/1000;

    H(1) = 0;
    H(2) = 11;
    H(3) = 25;
    H(4) = 47;
    H(5) = 53;
    H(6) = 79;
    H(7) = 90;
    H(8) = 100.001;
    H    = H*1000;  % convert to m

    g = 3.7159;
    C = -.019435;
    R = -g/C; %Condition Constant


    %%%%%%%%%%%%% DO NOT CHANGE BELOW THIS LINE %%%%%%%%%%%%%%%%%%%%%%%%%
    Tn = @(Tp,a,hn,hp) Tp+a*(hn-hp);

    Pn1 = @(Pp,Tn,Tp,C,hn,hp) Pp*(Tn/Tp)^(C/((Tn-Tp)/(hn-hp))); % Gradient
    Pn2 = @(Pp,Tn,Tp,C,hn,hp) Pp*exp(C/Tn*(hn-hp)) ;            % Isothermal

    Dn = @(Tn,Pn,R) Pn/(Tn*R);
    Hn = @(Hmax,Hin) (Hin>=Hmax)*Hmax+(Hin<Hmax)*Hin;

    L = 1; 

    while L<=length(H) && h>=H(L)
         T1 = Tn(T0,DT(L),Hn(H(L+1),h),H(L)); 
         if DT(L) == 0
            P0 = Pn2(P0,T1,T0,C,Hn(H(L+1),h),H(L));
         else
            P0 = Pn1(P0,T1,T0,C,Hn(H(L+1),h),H(L));  
         end
         T0 = T1;
         D0 = Dn(T0,P0,R);
         L  = L+1;    
end