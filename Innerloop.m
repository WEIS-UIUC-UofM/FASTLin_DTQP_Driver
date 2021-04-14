function [T1,U,X,Pitch_rate,F] = Innerloop(casefile,Pconstraint,WindFile)

clc;

mname = mfilename('fullpath');
[mpath,mname] = fileparts(mname);

% get linearized file path
savepath = fullfile(mpath,'SaveData',mname);

% get wind file
windpath = fullfile(mpath,'outb');

switch exist(savepath)
    case 0
        mkdir(savepath);
        disp(['Processed linear model location:',savepath]);
    case 7
        %disp(['Processed linear model location:',savepath]);
    otherwise
        error([savepath,' should be a directory.']);
end

save_name = strcat(casefile,'.mat');
LinModelFile = fullfile(savepath,save_name);

% load linearized files
if exist(LinModelFile,'file') == 2
    LinearModels = load(LinModelFile);
else
    % add openFAST matlab toolbox to path
    toolpath = strcat(mpath,'\matlab-toolbox-master');
    addpath(genpath(toolpath));
    FASTOutPath = strcat(mpath,'\SaveData\',casenames{ncase},'\linear\');
    ReduceModel = 0;
    Saveflag = 1;
    LinearModels = FASTLin_ProcessLinearModels(LinModelFile,FASTOutPath,ReduceModel,Saveflag,casenames{ncase});
end

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
    y_opsM(:,iSys) = LinearModels.SS_Ops(iSys).yop;
    
    % matrices
    Am(:,:,iSys) = LinearModels.P{iSys}.A;
    Bm(:,:,iSys) = LinearModels.P{iSys}.B;
    Cm(:,:,iSys) = LinearModels.P{iSys}.C;
    Dm(:,:,iSys) = LinearModels.P{iSys}.D;
    
end

% permute matrices
Am = permute(Am,[3,1,2]);
Bm = permute(Bm,[3,1,2]);
Cm = permute(Cm,[3,1,2]);
Dm = permute(Dm,[3,1,2]);

% interpolate system matrices based on wind speed
A_op = @(w) interp1(w_ops,Am,w,'pchip');
B_op = @(w) interp1(w_ops,Bm,w,'pchip');
C_op = @(w) interp1(w_ops,Cm,w,'pchip');
D_op = @(w) interp1(w_ops,Dm,w,'pchip');


%% load wind profile
if strcmp(WindFile,'072720_183300.mat') == 1
    WindFile = strcat(windpath,'\',WindFile);
    Wind_o = load(WindFile);
    Wind_speed = Wind_o.Chan.RtVAvgxh;
    tt = Wind_o.Chan.tt;
else
    WindFile = strcat(windpath,'\',WindFile);
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

%% disctretization points
opts.dt.nt = 2000;

%% W_fun
time = linspace(tt(1),tt(end),opts.dt.nt);
ppW = spline(tt,Wind_speed);
W_fun = @(t) ppval(ppW,t); % wind speed function

%% DxoDt_fun
ppDW = fnder(ppW);
DW_fun = @(t) ppval(ppDW,t);

Xo_pp = pchip(w_ops,x_opsM);
DXo_pp = fnder(Xo_pp);
DXo_fun = @(w) ppval(DXo_pp,w);

DXoDt_fun = @(t)((-DXo_fun(W_fun(t))).*DW_fun(t))';
%% Xo_fun
Xo_pp = pchip(w_ops,x_opsM);
Xo_fun = @(w) ppval(Xo_pp,w);

GS_pp = pchip(w_ops,x_opsM(5,:));
GS_fun = @(w)ppval(GS_pp,w);

GSn_pp = pchip(w_ops,-x_opsM(5,:));
GSn_fun = @(w)ppval(GSn_pp,w);

%% Uo_fun
Uo_pp = pchip(w_ops,u_opsM);
Uo_fun = @(w) ppval(Uo_pp,w);

GenTq_pp = pchip(w_ops,u_opsM(2,:));
GT_fun = @(w) ppval(GenTq_pp,w);

GenTqn_pp = pchip(w_ops,-u_opsM(2,:));
GTn_fun = @(w) ppval(GenTqn_pp,w) ;

BPitch_pp = pchip(w_ops,u_opsM(3,:));
BP_fun = @(w) ppval(BPitch_pp,w);
%% Yo_fun
Yo_pp = pchip(w_ops,y_opsM);
Yo_fun = @(w) ppval(Yo_pp,w);

%% yo(t)*uo(t)

GP = -u_opsM(2,:).*x_opsM(5,:);
GP_pp = pchip(w_ops,GP);
GP_fun = @(w) ppval(GP_pp,w);
%% state Bounds

% set upper and lower bounds for states
r = Xo_fun(w_ops);
indexat = @(expr, index) expr(index,:);

X0 = Xo_fun(W_fun(0));

GSr =max(r(5,:));
PTr = max(r(1,:));


ub = inf(size(max(r,[],2)));
ub(1) =   deg2rad(Pconstraint);
ub(5) = 0.7913+0.0001;
lb = -inf(size(max(r,[],2)));

UB1 = cell(101,1);
LB1 = cell(101,1);

for i = 1:101
    UB1{i} = @(t) ub(i) - indexat(Xo_fun(W_fun(t)),i);
    LB1{i} = @(t) lb(i) - indexat(Xo_fun(W_fun(t)),i);
end
X0(1) = ub(1);
X0(5) = ub(5);

%% DTQP options

opts.general.displevel = 2; % 0:none, 1:minimal, 2:verbose
opts.dt.defects = 'TR';     % ZO, EF, Huen, ModEF, TR, HS, RK4, PS
opts.dt.quadrature = 'CTR'; % CEF, CTR, CQHS, G, CC
opts.dt.mesh = 'ED';        % ED, LGL, CGL, USER
%opts.dt.meshr.method = 'RICHARDSON-DOUBLING';
opts.solver.function = 'quadprog';
opts.solver.tolerance = 1e-12;
opts.solver.maxiter = 5000;

%% Objective term

R1 = 1e-8;R2 = 1e+8;

% -(yl+yo)(ul+uo) + Ru^2
lx = 1;

% Ru^2
L(lx).left = 1;
L(lx).right = 1;
L(lx).matrix = diag([0,R1,R2]);
lx = lx+1;

% -yl*ul
L(lx).left = 1;
L(lx).right = 2;
Lmat = zeros(3,101); Lmat(2,5) = 1;
L(lx).matrix = -Lmat;
lx = lx+1;

% -ul*yo
L(lx).left = 0;
L(lx).right = 1;
L2mat = cell(1,3);L2mat{2} = @(t) GSn_fun(W_fun(t));
L(lx).matrix = L2mat;
lx = lx+1;

% -yl*uo
L3mat = cell(1,101); L3mat{5} = @(t) GTn_fun(W_fun(t));
L(lx).left = 0;
L(lx).right = 2;
L(lx).matrix = L3mat;
lx = lx+1;


TQmax = max(u_opsM(2,:));
GSmax = 0.7913;
SMat = zeros(101,1); SMat(5) = TQmax;

% -yo*uo
L(lx).left = 0;
L(lx).right = 0;
L(lx).matrix = {@(t) GP_fun(W_fun(t))};

% linear constraint on power
Z(1).linear(1).right = 1;
Z(1).linear(1).matrix = [0,GSmax,0]';
Z(1).linear(2).right = 2;
Z(1).linear(2).matrix = SMat;
Z(1).b = @(t)(15000 + 0.9655*TQmax*(GSmax + GSn_fun(W_fun(t))) + 0.9655*GSmax*(TQmax + GTn_fun(W_fun(t))));


%% Bounds
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

%% state bounds
ix = ix+1;
UB(ix).right = 2;
UB(ix).matrix = UB1;

LB(ix).right = 2;
LB(ix).matrix = LB1;


%% inital values
ix = ix+1;

load('X0_nominal.mat');

UB(ix).right = 4;
UB(ix).matrix =  X0_n-Xo_fun(W_fun(0));

LB(ix).right = 4;
LB(ix).matrix = X0_n-Xo_fun(W_fun(0));

%% setup

setup.A = TVmat2cell(A_op,time,W_fun,1);
setup.B = TVmat2cell(B_op,time,W_fun,1);
setup.d = TVmat2cell(DXoDt_fun,time,W_fun,0);
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
X = Xl+ Xo_fun(W_fun(T1))';
U = Ul + Uo_fun(W_fun(T1))';

Pitch = pchip(T1,rad2deg(U(:,3)));
Prate = fnder(Pitch);
Pr = @(t)ppval(Prate,t);
Pitch_rate = Pr(T1);
catch
   X = [];
   U = [];
   Pitch_rate = [];
end

%% Plot

plotflag = 1;

if plotflag
    figure
    
    subplot(3,2,1);
    plot(T1,U(:,1));
    title('Horizontal Wind Speed [m/s]');
    
    subplot(3,2,5);
    plot(T1,U(:,2)/1000);
    title('Generator Torque [KN-m]');
    ylim([1.8e+04 2e+04]);
    
    subplot(3,2,3);
    plot(T1,(U(:,3)));
    title('Blade Pitching Control [rad]');
    
    subplot(3,2,4);
    plot(T1,0.9655*X(:,5).*U(:,2)/1000);
    title('Gen speed * Gen Torque [kW]');
    
    subplot(3,2,6);
    plot(T1,X(:,5));
    title('Generator Speed [rad/s]');
    yticks([0.75 0.775 0.8]);
    
    subplot(3,2,2);
    tmp2 = plot(T1,rad2deg(X(:,1)));
    legend([tmp2],{'pitch'});
    title('Platform Tilt [deg]');
    
end

    
end



function A = TVmat2cell(f,time,W_fun,flag)

% function to convert nt*nx*nz matrix to nx*nx cell

% extract time 
t = time;


if flag
    w = W_fun(t);
else
    w = t;
end

% evaluate function
At = f(w);

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

