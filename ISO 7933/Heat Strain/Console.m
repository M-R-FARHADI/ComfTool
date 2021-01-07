clc
clear
%% Inputs from excel
dataset = xlsread('console.xlsx');
% weight set to default 75 kg
weight = dataset(1,1);
% height set to default 1.8 m
height = dataset(2,1);
% DRINK set to default 1
DRINK = dataset(6,1);
% Adu calculated by default with height and weight
Adu = 0.202 * (weight ^ 0.425) * (height ^ 0.725);
% spHeat (=c_sp) calculated by default with Adu and weight
spHeat = 57.83 * weight / Adu;
% Swp set to default 0
SWp = 0;
% Tcr set to default 36.8
Tcr = dataset(3,1);
% Tre set to Tcr
Tre = Tcr;
% Tcreq set equal Tcr
Tcreq = Tcr;
% Tsk set to default 34.1
Tsk = dataset(4,1);
% The user must make sure that, at this point in the programme,
% the following parameters are available.
% Standard values must be replaced by actual values if necessary.
% The water replacement is supposed to be sufficient so that the
% workers can drink freely (DRINK=1), otherwise the value DRINK=0 must be used
% following lines from code in standard
SWtot = 0;
TskTcrwg = .3;
Dmax50 = 0.075 * weight * 1000;
Dmax95 = 0.05 * weight * 1000;

% EXPONENTIAL AVERAGING CONSTANTS
% Core temperature as a function of
% the metabolic rate: time constant: 10 minutes
ConstTeq = exp(-1 / 10);
% Skin Temperature: time constant: 3 minutes
ConstTsk = exp(-1 / 3);
% Sweat rate: time constant: 10 minutes
ConstSW = exp(-1 / 10);
% Duration set to standard of 480 minutes
% the duration of the work sequence in minutes
Duration = 480;

% setting default values originally in the for loop
% air temperature degrees celsius
Ta = dataset(5,1);
% Tr set by default to Ta
Tr = Ta;
Pa = dataset(8,1);
%HR = 50; % normal barometric pressure in [Pa]
%Pa = HR/1000/(0.622+HR/1000)*pb/1000;
% air velocity metres per second
Va = dataset(7,1);
% metabolic rate Watts per square meter
Met = dataset(9,1);
% effective mechanical power Watts per square metre
Work = dataset(10,1);
% posture set by default to 2 = standing
%posture = 1 sitting, = 2 standing, = 3 crouching
posture = dataset(11,1);
% Icl set by default to 0.5 clo
% static thermal insulation clo
Icl = dataset(12,1);
% imst set by default to 0.38
%static moisture permeability index dimensionless
imst = 0.38;
% Ap set by default to 0.54
% fraction of the body surface covered by the
% reflective clothing dimensionless
Ap = 0.54;
% emissivity of the reflective clothing dimensionless (by default: Fr=0.97)
% Fr set by default to 0.97
Fr = 0.97;
% accl set by default to 100 = acclimated subject
%code =100 if acclimatised subject, 0 otherwise
accl = 100;
THETA = dataset(13,1); %angle between walking direction and wind direction degrees
defdir = 0; %code =1 if walking direction entered, 0 otherwise
if THETA ~= 0
    defdir = 1;
end
Walksp = dataset(13,1);  %walking speed metres per second
if (Walksp == 0)    
    defspeed = 0;  
else 
    defspeed = 1;
end
[Tre,SWtotg,Dlimtre,Dlimloss50,Dlimloss95,...
    Cres,Eres,Ep,SWp,Texp,Tskeq,Tsk,wp] = iso7933(accl, posture,...
    Ta, Pa, Tr, Va, Met, Icl, THETA, Walksp, Duration, weight,...
    DRINK, Adu, spHeat, SWp, Tre, Tcr, Tsk, Tcreq, Work,...
    imst, Ap, Fr, defspeed, defdir)