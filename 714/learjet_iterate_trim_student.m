%==============================================================================
%   M-File to set up parameters on learjet longitudinal trim simulation
%==============================================================================
%       Data originates from Roskam "Airplane Flight Dynamics" book, 
%       Airplane G Appendix B Page 522-523
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  input below the solution from solving approx trim problem eqn 4.45 page 199
clear; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
g = 32.2;

%  Geometric Data
Sref = 232;
cbar = 7.04;
bref = 34;
U1 = 170.0;  % ft/sec
it = 0;  % thrust angle
dT=0;  % thrust offset
cosit = cos(it);
sinit = sin(it);
qbar = 34.3;

%Mass Properties

Wt = 13000;
mass = Wt/g;
Ixxb = 28000;
Iyyb = 18800;
Izzb = 47000;
Ixzb = 1300;
inertia_matrix = [Ixxb 0 -Ixzb; 0 Iyyb 0;Ixzb 0 Izzb];

%   Aero data Airplane G approach condition
CL0 = 1.2;
CLA = 5.04;
CLADOTHAT = 1.6;
CLQHAT = 4.1;
CLDE = .40;
CLDETAB = CLDE/6.;

CM0 = .047;
CMA = -.66   ;
CMADOTHAT = -5.0;
CMQHAT = -13.5;
CMDE = -0.98;
CMDETAB = CMDE/6;
dPMdT = 0;
%  drag polar parameters below
CD0BAR = .01867;
CDK = .18824;
%Linear Drag model parameters below (can choose between 2 models in drag block)
CD0=.0431; %(.1635 modified from value in Roskam of .0431 to match CD1 from Roskam)
CDA=1.06;
CDDE = 0;
CDDETAB = 0;
%
%  set up iterative trim solution
%
% Matrix problem is AMAT*[ alphatrim, detrim ] = [ CL1-CL0-CT1*sin(it+alphatrim), -CM0+CT1*dT/cbar ]
CL1 = Wt/(qbar*Sref);
AMAT = [CLA, CLDE ;
        CMA, CMDE];
iAMAT = inv(AMAT);
rhs = @(alpha, T1) [ CL1-CL0-(T1/(qbar*Sref))*sin(it+alpha); -CM0+(T1/(qbar*Sref))*dT/cbar ] ;
rhsf = [ CL1-CL0; -CM0 ] ;
thrust = @(alpha, de) ((CD0+CDA*alpha+CDDE*de)*qbar*Sref)/(cos(it+alpha));
varhat = iAMAT*rhsf;
alphahat= varhat(1);
dehat = varhat(2);
thrusthat = thrust(alphahat,dehat);
err = 1;
tol = 1e-4;
while err > tol
    varnew = iAMAT*rhs(alphahat,thrusthat);
    alphanew = varnew(1) ;
    denew = varnew(2);
    thrustnew = thrust(alphahat,dehat);
%     err = abs((alphanew - alphahat)/alphanew) + abs((denew-dehat)/denew) + abs((thrustnew-thrusthat)/thrustnew)
    err = abs((alphanew - alphahat));
    alphahat = alphanew;
    dehat = denew;
    thrusthat = thrustnew;
end
thrust = thrusthat;
alphatrim = alphahat*(180/pi);
detrim = dehat*(180/pi);

%   Set up IC's on Euler Angle block.
alphatrim=alphatrim*pi/180; %need radians for sin and cos here
alpha=alphatrim; %make sure you have converted to radians for sin and cos before here
hic = 0;
Uic = U1*cos(alpha);
Vic = 0;
Wic = U1*sin(alpha);
thetaic = alpha;

detab = 0;
psiic = 0;
phiic = 0;
%%%%%%      Uic is HAND FIXED
Vtrue_ic=sqrt(Uic^2+Vic^2+Wic^2);
%%%%%%

%   Run the simulation.  
%       Length of run-time is set by the value of TIMESPAN
sim ('learjet')

%   Plot the output
%plot6