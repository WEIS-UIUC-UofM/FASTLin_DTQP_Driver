% Innerloop.m
% Solves the inner-loop optimal control problem using direct transcription
% and quadratic programming through the DTQP tool.
% Primarily this file sets up the problem elements for DTQP
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------

function [T,U,X,Pitch_rate,F] = Innerloop(LinearModels,Pconstraint,WindFile,PlotFlag,opts)

% Set length
nl = length(LinearModels.P);

% go through each linear model
for iSys = 1:nl
    
    % extract
    sys = LinearModels.P{iSys};
    
    % number of states
    nStates(iSys) = size(LinearModels.P{iSys}.A,1);
    
    % wind speed
    w_ops(iSys,1) = LinearModels.WindSpeed(iSys);
    
    % operating points
    x_opsM(:,iSys) = LinearModels.SS_Ops(iSys).xop;
    u_opsM(:,iSys) = LinearModels.SS_Ops(iSys).uop;
    
    % matrices
    Am(:,:,iSys) = LinearModels.P{iSys}.A;
    Bm(:,:,iSys) = LinearModels.P{iSys}.B;
    
end

% permute matrices
Am = permute(Am,[3,1,2]);
Bm = permute(Bm,[3,1,2]);

% interpolate system matrices based on wind speed
A_op = @(w) interp1(w_ops,Am,w,'pchip');
B_op = @(w) interp1(w_ops,Bm,w,'pchip');

%% load wind profile
if strcmp(WindFile,'072720_183300.mat') == 1

    Wind_o = load(WindFile);
    Wind_speed = Wind_o.Chan.RtVAvgxh;
    tt = Wind_o.Chan.tt;
else

    Wind_o = load(WindFile);
    idx = find(contains(Wind_o.ChanName,'RtVAvgxh'));
    Wind_speed = Wind_o.Channels(:,idx);
    tt = Wind_o.Channels(:,1);    
end

filterflag = 1;

% filter
if filterflag
    t_f     = 1;
    dt      = tt(2) - tt(1);
    
    nb      = floor(t_f/dt);
    b       = ones(nb,1)/nb;
    
    Wind_speed     = filtfilt(b,1,Wind_speed);
end

%% W_fun
% interpolate the discrete values of the wind speed into a continuous
% function
time = linspace(tt(1),tt(end),opts.dt.nt)';
ppW = spline(tt,Wind_speed);
W_fun = @(t) ppval(ppW,t); % wind speed function

%% DxoDt_fun
% LPV model correction term for time varying nature of the operating point
% (Dxo/Dt)

% extract details from piece-wise polynomial by breaking it apart
[breaks,coefs,l,k,d] = unmkpp(ppW);

% make the polynomial that describes the derivative
ppDW = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
DW_fun = @(t) ppval(ppDW,t);

Xo_pp = pchip(w_ops,x_opsM);

% extract details from piece-wise polynomial by breaking it apart
[breaks,coefs,l,k,d] = unmkpp(Xo_pp);

% make the polynomial that describes the derivative
DXo_pp = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
DXo_fun = @(w) ppval(DXo_pp,w);

DXoDt_fun = @(t)((-DXo_fun(W_fun(t'))).*DW_fun(t'))';

%% Interpolate the discrete operating point values into continuous functions
Disc2Cont
%% state Bounds

% set upper and lower bounds for states
r = Xo_fun(w_ops);

% function to find the values of anaonymous functions at specific index
indexat = @(expr, index) expr(index,:); 

% set upper bound values
ub = inf(size(max(r,[],2)));

% PtfmPitch constraint
ub(1) =   deg2rad(Pconstraint);

% GenSpeed constraint
ub(5) = 0.7913+0.0001; 

% lower bound values
lb = -inf(size(max(r,[],2)));

% initialize UB,LB cell arrays
UB1 = cell(101,1);
LB1 = cell(101,1);

% 
for i = 1:101
    UB1{i} = @(t) ub(i) - indexat(Xo_fun(W_fun(t)),i);
    LB1{i} = @(t) lb(i) - indexat(Xo_fun(W_fun(t)),i);
end


%% DTQP options

opts.dt.defects = 'TR';     % ZO, EF, Huen, ModEF, TR, HS, RK4, PS
opts.dt.quadrature = 'CTR'; % CEF, CTR, CQHS, G, CC
opts.dt.mesh = 'ED';        % ED, LGL, CGL, USER
%opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
opts.solver.function = 'quadprog';
opts.solver.tolerance = 1e-12;
opts.solver.maxiter = 5000;

%% Objective term
Lagrange

%% Control Bounds
ix = 1;

% remove the offset for control bounds

UB(ix).right = 1;
UB(ix).matrix = {@(t) W_fun(t)-W_fun(t);
    @(t) max(u_opsM(2,:))-GT_fun(W_fun(t));
    @(t) max(u_opsM(3,:))-BP_fun(W_fun(t))};

LB(ix).right = 1;
LB(ix).matrix = {@(t) W_fun(t)-W_fun(t);
    @(t) min(u_opsM(2,:))-GT_fun(W_fun(t));
    @(t) min(u_opsM(3,:))-BP_fun(W_fun(t))};

%% State Bounds
ix = ix+1;

UB(ix).right = 2;
UB(ix).matrix = UB1;

LB(ix).right = 2;
LB(ix).matrix = LB1;


%% Inital Values
ix = ix+1;

load('X0_nominal.mat');

UB(ix).right = 4;
UB(ix).matrix =  X0_n-Xo_fun(W_fun(0));

LB(ix).right = 4;
LB(ix).matrix = X0_n-Xo_fun(W_fun(0));

%% DTQP setup

setup.A = TVmat2cell(@(t)A_op(W_fun(t)),time);
setup.B = TVmat2cell(@(t)B_op(W_fun(t)),time);
setup.d = TVmat2cell(DXoDt_fun,time);
setup.L = L;
setup.Z = Z;
setup.UB = UB;
setup.LB = LB;
setup.t0 = 0;
setup.tf = time(end);

%% solve
[T1,Ul,Xl,P,F,in,opts] = DTQP_solve(setup,opts);

%% Add the offset values for controls and states
try
T = T1;    
X = Xl+ Xo_fun(W_fun(T1))';
U = Ul + Uo_fun(W_fun(T1))';

% evaluate pitch rate
Pitch = pchip(T1,rad2deg(U(:,3)));
[breaks,coefs,l,k,d] = unmkpp(Pitch);
Prate = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
Pitch_rate = ppval(Prate,T1);


if PlotFlag
    % plot
    IL_plot
end

catch
   X = [];
   U = [];
   Pitch_rate = [];
end

%% Plot
    
end



function A = TVmat2cell(f,time)

% function to convert nt*nx*nz matrix to nx*nx cell

% evaluate function
At = f(time);

% get sizes
m = size(At,2);
n = size(At,3);

% initialize cell
A = cell(m,n);

for i = 1:m
    for j = 1:n
        A{i,j} = @(t) interp1(time,At(:,i,j),t);
    end
end

end


