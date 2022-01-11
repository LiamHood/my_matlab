%==============================================================================
%   M-File to set up parameters on learjet longitudinal trim simulation
%==============================================================================
%       Data originates from Roskam "Airplane Flight Dynamics" book, 
%       Airplane G Appendix B Page 522-523
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  input below the solution from solving approx trim problem eqn 4.45 page 199
alphatrim = 4        %  input initial trim alpha in degrees
detrim =    0       %  input initial trim elevator in degrees 
thrust =    1000     % input initial trim thrust in lbf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
g = 32.2

%  Geometric Data
Sref = 232
cbar = 7.04
bref = 34.
U1 = 170.0  % ft/sec
it = 0  % thrust angle
cosit = cos(it)
sinit = sin(it)

%Mass Properties

Wt = 13000.
mass = Wt/g
Ixxb = 28000
Iyyb = 18800
Izzb = 47000
Ixzb = 1300
inertia_matrix = [Ixxb 0 -Ixzb; 0 Iyyb 0;Ixzb 0 Izzb]

%   Aero data Airplane G approach condition
CL0 = 1.2
CLA = 5.04
CLADOTHAT = 1.6
CLQHAT = 4.1
CLDE = .40
CLDETAB = CLDE/6.

CM0 = .047
CMA = -.66   
CMADOTHAT = -5.0
CMQHAT = -13.5
CMDE = -0.98
CMDETAB = CMDE/6
dPMdT = 0
%  drag polar parameters below
CD0BAR = .01867
CDK = .18824
%Linear Drag model parameters below (can choose between 2 models in drag block)
CD0=.1635 %(modified from value in Roskam of .0431 to match CD1 from Roskam)
CDA=1.06
CDDE = 0
CDDETAB = 0

%   Set up IC's on Euler Angle block.
alpha = alphatrim*pi/180.0  %  convert from deg to radians
hic = 0
Uic = U1*cos(alpha)
Vic = 0
Wic = U1*sin(alpha)
thetaic = alpha

detab = 0.
psiic = 0
phiic = 0
%%%%%%      Uic is HAND FIXED
Vtrue_ic=sqrt(Uic^2+Vic^2+Wic^2)
%%%%%%

%   Run the simulation.  
%       Length of run-time is set by the value of TIMESPAN
sim ('learjet')

%   Plot the output
%plot6