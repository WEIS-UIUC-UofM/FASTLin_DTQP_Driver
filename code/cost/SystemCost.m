function [Mturbine,Cturbine] = SystemCost(Mfactor)
%% System plant parameters

% Rotor diameter
Drot = 117*2;

% Tower height
Lhub = 150;

% Turbine rating
Ptur = 15;

% Number of blades
nblade = 3;

% Diameter of chain
d = 185/1000;

%  Number of mooring lines
nmoor = 3;

% Length of mooring lines
Lmoor = 850;

%% Blade cost

km = 0.5;b = 2.47;kc = 14.6;

Mblade = km*(0.5*Drot)^b;
Cblade = kc*Mblade;

%% Generator cost

km = 2300;b = 3400;kc=12.4;

Mgen = km*Ptur + b;
Cgen = kc*Mgen;

%% Tower cost

km = 19.828;b = 2.0282;kc=2.9;

Mtower = km*(Lhub)^b;
Ctower = kc*Mtower;

%% Hub systems

% Hub
km = 2.3; b = 1320; kc = 3.9;

Mhub = km*Mblade + b;
Chub = kc*Mhub;

% Pitch system
km = 0.1295; b1 = 491.31; b2 = 555; h = 0.328; kc = 22.1;

Mbearing = nblade*km*Mblade+b1;
Mpitch = Mbearing*(1+h) + b2;
Cpitch = kc*Mpitch;

% Spinner
km = 15.5; b = -980; kc = 11.1;

Mspin = km*Drot + b;
Cspin = kc*Mspin;

%% Nacelle

% Low speed shaft
km = 13; b1 = 0.65; b2 = 775; kc = 11.9;

Mls = km*(Mblade*Ptur)^b1 + b2;
Cls = kc*Mls;

% High speed shaft
km = 198.94; kc = 6.8;

Mhs = km*Ptur;
Chs = kc*Mhs;

% Yaw system
km = 0.00135; b = 3.314; kc = 8.3;

Myaw = km*(Drot)^b;
Cyaw = kc*Myaw;

% Hydraulic cooling
km = 80; kc = 124;

Mhc = km*Ptur;
Chc = kc*Mhc;

% Transformer
km = 1915; b = 1910; kc = 18.8;

Mtrans = km*Ptur + b;
Ctrans = kc*Mtrans;

% Cabling
kc = 41850;

Ccab = kc*Ptur;

% Control System
kc = 21150;

Ccont = kc*Ptur;

% Nacelle cover
km = 1.2817; b = 428.19; kc = 5.7;

Mcover = km*Ptur + b;
Ccover = kc*Mcover;

%% Mooring

Mmoor = nmoor*Lmoor*19.9*1000*d^2;
Cmoor = 3.415*10^4*d^2;

%% Platform

if varargin < 1
    Mfactor = 0.2:0.1:1.2;
end

% component weights
M_steel = 3914*1000;
M_ballast = 2540*1000;
M_outfit = 100*1000;
M_sw = 11300*1000+16000;

% cost per tonne of components
C_steel = 3000*1.22; % euros to dollars
C_ballast = 150;
C_outfit = 7250;

Mptfm = M_steel + M_ballast + M_outfit + M_sw;
Cptfm = (M_steel/1000*C_steel + M_ballast/1000*C_ballast + M_outfit/1000*C_outfit )*Mfactor;

%% Aggregate cost

% Hub system
Mhubsys = Mhub + Mpitch + Mspin;
Chubsys = Chub + Cpitch + Cspin;

% Rotor system
Mrotor = nblade*Mblade + Mhubsys;
Crotor = nblade*Cblade + Chubsys;

% Nacelle system
Mnacelle = Mls + Mhs + Mgen + Myaw + Mhc + Mtrans + Mcover;
Cnacelle = Cls + Chs + Cgen + Cyaw + Chc + Ctrans + Ccab + Ccont + Ccover;

% Turbine
Mturbine = Mrotor + Mnacelle + Mtower + Mmoor + Mptfm;
Cturbine = Crotor + Cnacelle + Ctower + Cmoor + Cptfm;

%% Simplified cost

Cturbine_simple = (nblade*Cblade + Ctower + Cgen + Cmoor);


%% Old costs

C_tower = (0.3973*pi*(Drot/2)^2*Lhub-1414)*1.5;
C_blades = (0.4019*(Drot/2)^3-955.24)*3;
C_nacelle = (11.537*Ptur*1000+3849.7);
C_gen = Ptur*1000*219.33;
C_mooring = 0.42*Lmoor*2437;

C_turbine_old = C_tower + C_blades + C_nacelle + C_gen + C_mooring;

end