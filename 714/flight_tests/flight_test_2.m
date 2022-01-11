close all; clear; clc;

%% Tests
load("Data_ft2.mat")
% filters = Data747approachlight.VVI__fpm<500 & Data747approachlight.VVI__fpm>-500 ...
%     & Data747approachlight.Vtrue_ktas < 132 & Data747approachlight.Vtrue_ktas > 128 ...
%     & Data747approachlight.alt_1ftmsl < 5200 & Data747approachlight.alt_1ftmsl > 4800;
% velocity = Data747approachlight.Vtrue_ktas(filters);
% figure
% plot(velocity,'.')
xJa = find(Data747approachlight.Vtrue_ktas <= 129.032 & Data747approachlight.Vtrue_ktas >= 129.030,1);
xJaalpha = Data747approachlight.alpha__deg(xJa);
xJade = Data747approachlight.elev1__deg(xJa);
xJaih = 0;

% filters = Data747cruise.VVI__fpm<500 & Data747cruise.VVI__fpm>-500 ...
%     & Data747cruise.Vtrue_ktas < 400.4 & Data747cruise.Vtrue_ktas > 396.4...
%     & Data747cruise.alt_1ftmsl < 20200 & Data747cruise.alt_1ftmsl > 19800;
% velocity = Data747cruise.Vtrue_ktas(filters);
% figure
% plot(velocity,'.')
xJc = find(Data747cruise.Vtrue_ktas <= 398.785 & Data747cruise.Vtrue_ktas >= 398.783,1);
xJcalpha = Data747cruise.alpha__deg(xJc);
xJcde = Data747cruise.elev1__deg(xJc);
xJcih = 0;

% filters = DataCessna.VVI__fpm<500 & DataCessna.VVI__fpm>-500 ...
%     & DataCessna.Vtrue_ktas < 65.4032 & DataCessna.Vtrue_ktas > 61.4032...
%     & DataCessna.alt_1ftmsl < 5200 & DataCessna.alt_1ftmsl > 4800;
% velocity = DataCessna.Vtrue_ktas(filters);
% figure
% plot(velocity,'.')
xAa = find(DataCessna.Vtrue_ktas <= 65.040 & DataCessna.Vtrue_ktas >= 65.0390,1);
xAaalpha = DataCessna.alpha__deg(xAa);
xAade = DataCessna.elev1__deg(xAa);
xAaih = 0;

% filters = DataCessna.VVI__fpm<500 & DataCessna.VVI__fpm>-500 ...
%     & DataCessna.Vtrue_ktas < 132.2992 & DataCessna.Vtrue_ktas > 128.2992 ...
%     & DataCessna.alt_1ftmsl < 5200 & DataCessna.alt_1ftmsl > 4800;
% velocity = DataCessna.Vtrue_ktas(filters);
% figure
% plot(velocity,'.')
xAc = find(DataCessna.Vtrue_ktas <= 129.927 & DataCessna.Vtrue_ktas >= 129.926,1);
xAcalpha = DataCessna.alpha__deg(xAc);
xAcde = DataCessna.elev1__deg(xAc);
xAcih = 0;

% filters = DataF4.VVI__fpm<500 & DataF4.VVI__fpm>-500 ...
%     & DataF4.Vtrue_ktas < 138.1600 & DataF4.Vtrue_ktas > 134.1600...
%     & DataF4.alt_1ftmsl < 5200 & DataF4.alt_1ftmsl > 4800;
% velocity = DataF4.Vtrue_ktas(filters);
% figure
% plot(velocity,'.')
xIa = find(DataF4.Vtrue_ktas <= 138.13 & DataF4.Vtrue_ktas >= 138.11,1);
xIaalpha = DataF4.alpha__deg(xIa);
xIade = DataF4.elev1__deg(xIa);
xIaih = 0;

% filters = DataF4.VVI__fpm<500 & DataF4.VVI__fpm>-500 ...
%     & DataF4.Vtrue_ktas < 520.5920 & DataF4.Vtrue_ktas > 516.5920 ...
%     & DataF4.alt_1ftmsl < 35200 & DataF4.alt_1ftmsl > 34800;
% velocity = DataF4.Vtrue_ktas(filters);
% figure
% plot(velocity,'.')
xIc = find(DataF4.Vtrue_ktas <= 517.246 & DataF4.Vtrue_ktas >= 517.244,1);
xIcalpha = DataF4.alpha__deg(xIc);
xIcde = DataF4.elev1__deg(xIc);
xIcih = 0;


%% TAS to IAS
% Aa, Ac, Ja, Jc, Ia, Ic
tas = [107.1, 220.1, 221, 673, 230, 876];
rhos = 23.77e-4;
rho1 = [20.48e-4,20.48e-4,20.48e-4,20.48e-4,20.48e-4,20.48e-4];
ps = 14.696;
p1 = [12.228,12.228,12.228,12.228,12.228,12.228];
as = 1115.64;
a1 = [1097.08,1097.08,1097.08,1097.08,1097.08,1097.08];
for ii = 1:length(tas)
    ias(ii) = TAStoIAS(tas(ii),rhos,rho1(ii),ps,p1(ii),as,a1(ii));
end
disp(tas)
disp(ias)
disp(tas*.592)
disp(ias*.592)

%% Cessna 172 Appendix A
Aac0.D = .0605;
Acc0.D = .0270;
Aac0.L = .807;
Acc0.L = .307;
Aac0.m = .09;
Acc0.m = .04;

Aaca.D = .547;
Acca.D = .121;
Aaca.L = 4.41;
Acca.L = 4.41;
Aaca.m = -.611;
Acca.m = -.613;

Aacih.D = 0;
Accih.D = 0;
Aacih.L = 0;
Accih.L = 0;
Aacih.m = 0;
Accih.m = 0;

Aacde.D = 0;
Accde.D = 0;
Aacde.L = .43;
Accde.L = .43;
Aacde.m = -1.029;
Accde.m = -1.122;

Aaalpha = 4;
Acalpha = 0;
Aaih = 0;
Acih = 0;
Aade = 4;
Acde = 0;
[AaC.D, AaC.L, AaC.m] = cL_linear(Aac0, Aaca, Aacih, Aacde, Aaalpha*pi/180, Aaih*pi/180, Aade*pi/180);
[AcC.D, AcC.L, AcC.m] = cL_linear(Acc0, Acca, Accih, Accde, Acalpha*pi/180, Acih*pi/180, Acde*pi/180);
[AaCx.D, AaCx.L, AaCx.m] = cL_linear(Aac0, Aaca, Aacih, Aacde, xAaalpha*pi/180, xAaih*pi/180, xAade*pi/180);
[AcCx.D, AcCx.L, AcCx.m] = cL_linear(Acc0, Acca, Accih, Accde, xAcalpha*pi/180, xAcih*pi/180, xAcde*pi/180);
disp("Cessna 172 Appendix A")
disp("Approach Book")
disp(AaC)
disp("Approach Test")
output(DataCessna, xAa)
disp(AaCx)
disp("Cruise Book")
disp(AcC)
disp("Cruise Test")
output(DataCessna, xAc)
disp(AcCx)

%% 747-400 Appendix J
Jac0.D = .0751;
Jcc0.D = .0164;
Jac0.L = .92;
Jcc0.L = .21;
Jac0.m = 0;
Jcc0.m = 0;

Jaca.D = 1.13;
Jcca.D = .2;
Jaca.L = 5.67;
Jcca.L = 4.4;
Jaca.m = -1.45;
Jcca.m = -1.00;

Jacih.D = 0;
Jccih.D = 0;
Jacih.L = .75;
Jccih.L = .7;
Jacih.m = -3.0;
Jccih.m = -2.7;

Jacde.D = 0;
Jccde.D = 0;
Jacde.L = .36;
Jccde.L = .32;
Jacde.m = -1.4;
Jccde.m = -1.3;

Jaalpha = 8.5;
Jcalpha = 2.5;
Jaih = 0;
Jcih = 0;
Jade = 4;
Jcde = 0;
[JaC.D, JaC.L, JaC.m] = cL_linear(Jac0, Jaca, Jacih, Jacde, Jaalpha*pi/180, Jaih*pi/180, Jade*pi/180);
[JcC.D, JcC.L, JcC.m] = cL_linear(Jcc0, Jcca, Jccih, Jccde, Jcalpha*pi/180, Jcih*pi/180, Jcde*pi/180);
[JaCx.D, JaCx.L, JaCx.m] = cL_linear(Jac0, Jaca, Jacih, Jacde, xJaalpha*pi/180, xJaih*pi/180, xJade*pi/180);
[JcCx.D, JcCx.L, JcCx.m] = cL_linear(Jcc0, Jcca, Jccih, Jccde, xJcalpha*pi/180, xJcih*pi/180, xJcde*pi/180);
disp("747-400 Appendix J")
disp("Approach Book")
disp(JaC)
disp("Approach Test")
output(Data747approachlight, xJa)
disp(JaCx)
disp("Cruise Book")
disp(JcC)
disp("Cruise Test")
output(Data747cruise, xJc)
disp(JcCx)

%% F-4 Appendix I
Iac0.D = .0269;
Icc0.D = .0205;
Iac0.L = .43;
Icc0.L = .1;
Iac0.m = .02;
Icc0.m = .025;

Iaca.D = .555;
Icca.D = .3;
Iaca.L = 2.8;
Icca.L = 3.75;
Iaca.m = -.098;
Icca.m = -.4;

Iacih.D = -.14;
Iccih.D = -.10;
Iacih.L = .24;
Iccih.L = .4;
Iacih.m = -.322;
Iccih.m = -.580;

Iacde.D = 0;
Iccde.D = 0;
Iacde.L = 0;
Iccde.L = 0;
Iacde.m = 0;
Iccde.m = 0;

Iaalpha = 11.7;
Icalpha = 2.6;
Iaih = 0;
Icih = 0;
Iade = 0;
Icde = 0;
[IaC.D, IaC.L, IaC.m] = cL_linear(Iac0, Iaca, Iacih, Iacde, Iaalpha*pi/180, Iaih*pi/180, Iade*pi/180);
[IcC.D, IcC.L, IcC.m] = cL_linear(Icc0, Icca, Iccih, Iccde, Icalpha*pi/180, Icih*pi/180, Icde*pi/180);

[IaCx.D, IaCx.L, IaCx.m] = cL_linear(Iac0, Iaca, Iacih, Iacde, xIaalpha*pi/180, xIaih*pi/180, xIade*pi/180);
[IcCx.D, IcCx.L, IcCx.m] = cL_linear(Icc0, Icca, Iccih, Iccde, xIcalpha*pi/180, xIcih*pi/180, xIcde*pi/180);
disp("F-4 Appendix I")
disp("Approach Book")
disp(IaC)
disp("Approach Test")
output(DataF4, xIa)
disp(IaCx)
disp("Cruise Book")
disp(IcC)
disp("Cruise Test")
output(DataF4, xIc)
disp(IcCx)

%% functions
function [ias] = TAStoIAS(tas,rhos,rho1,ps,p1,as,a1)
    M = tas/a1 ;
    if M < .3
        ias = tas/sqrt(rhos/rho1);
    else
        gamma = 1.4;
        p0 = p1*(1 + (gamma-1)/2*M^2)^(gamma/(gamma-1));
        ias = sqrt( (2*as^2)/(gamma-1)*( (((p0-p1)/ps)+1)^((gamma-1)/gamma)-1));
    end
end
 
function [CD, CL, Cm] = cL_linear(c0, ca, cih, cde, alpha, ih, de)
    coefficients = [c0.D, ca.D, cih.D, cde.D;
                    c0.L, ca.L, cih.L, cde.L;
                    c0.m, ca.m, cih.m, cde.m];
    state = [1; alpha; ih; de];
    C = coefficients*state;
    CD = C(1);
    CL = C(2);
    Cm = C(3);
end

function output(data, x)
    disp(["start time",max(data.real_time(data.real_time<(data.real_time(x)-5)))])
    disp(["end time",data.real_time(x)])
    disp(["ias",data.Vind_kias(x)])
    disp(["tas",data.Vtrue_ktas(x)])
    disp(["alt",data.alt_1ftmsl(x)])
    disp(["AoA",data.alpha__deg(x)])
    disp(["weight",data.curnt___lb(x)])
    disp(["elevator", data.elev1__deg(x)])
    disp(["M",data.M_ftlb(x)])
    disp(["CL",data.cltotal(x)])
    disp(["CD",data.cdtotal(x)])
end
