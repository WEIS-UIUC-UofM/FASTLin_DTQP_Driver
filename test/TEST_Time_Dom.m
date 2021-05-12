%--------------------------------------------------------------------------
% TEST_Time_Dom.m
% Time domain comparison between LPV model and NL openfast
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
%close all; 
clear; clc;close all;

% plot flag
debugPlotFlag = false;

% reduce the number of time points for defining the inputs (N below)
reducedInputFlag = false;

% load linear models
LinearModels = load('pd_1.0_linear.mat');

% extract
P = LinearModels.P;

% number of linear models
nl = length(P);

% size of the matrices
sizeA = size(P{1}.A);
sizeB = size(P{1}.B);
sizeC = size(P{1}.C);
sizeD = size(P{1}.D);

% initialize sparsity pattern matrices
Asparsity = zeros(sizeA);
Bsparsity = zeros(sizeB);
Csparsity = zeros(sizeC);
Dsparsity = zeros(sizeD);

% initialize storage matrices for each wind speed
Aw = zeros([nl size(Asparsity)]);
Bw = zeros([nl size(Bsparsity)]);
Cw = zeros([nl size(Csparsity)]);
Dw = zeros([nl size(Dsparsity)]);

% go through each linear model
for iSys = 1:nl

    % extract
    Psys = P{iSys};
    SS_Ops_sys = LinearModels.SS_Ops(iSys);

    % state indices
    StateNames = Psys.StateName;
    iGenSpeed = find(strcmpi('ED First time derivative of Variable speed generator DOF (internal DOF index = DOF_GeAz), rad/s',StateNames)); % generator speed
    iPltPitch = find(strcmpi('ED Platform pitch tilt rotation DOF (internal DOF index = DOF_P), rad',StateNames)); % platform pitch

    % input indices
    InputNames = Psys.InputName;
    iWindSpeed = find(strcmpi('IfW WindSpeedHor',InputNames)); % wind speed
    iGenTorq = find(strcmpi('ED GenTorq',InputNames)); % generator torque
    iBldPitch = find(strcmpi('ED BldPitchCommand',InputNames)); % blade pitch
    indOuts = [iWindSpeed,iGenTorq,iBldPitch];
%     indOuts = [iGenTorq,iWindSpeed,iBldPitch];

    % output indices
    OutputNames = Psys.OutputName;

    % number of states
    nStates(iSys) = size(Psys.A,1);

    % wind speed
    w_ops(iSys,1) = SS_Ops_sys.uop(iWindSpeed);

    % operating points
    x_opsM(:,iSys) = SS_Ops_sys.xop;
    u_opsM(:,iSys) = SS_Ops_sys.uop(indOuts);
    y_opsM(:,iSys) = SS_Ops_sys.yop;

    % matrices
    Aw(iSys,:,:) = Psys.A;
    Bw(iSys,:,:) = Psys.B(:,indOuts);
    Cw(iSys,:,:) = Psys.C;
    Dw(iSys,:,:) = Psys.D;

    % sparsity patterns
    Asparsity = Asparsity|shiftdim(Aw(iSys,:,:));
    Bsparsity = Bsparsity|shiftdim(Bw(iSys,:,:));
    Csparsity = Csparsity|shiftdim(Cw(iSys,:,:));
    Dsparsity = Dsparsity|shiftdim(Dw(iSys,:,:));

end

% flatten sparsity patterns and nonzero data
logicalA = Asparsity(:);
Aw_v = reshape(Aw,iSys,[]);
Aw_v = Aw_v(:,logicalA);

logicalB = Bsparsity(:);
Bw_v = reshape(Bw,iSys,[]);
Bw_v = Bw_v(:,logicalB);

logicalC = Csparsity(:);
Cw_v = reshape(Cw,iSys,[]);
Cw_v = Cw_v(:,logicalC);

logicalD = Dsparsity(:);
Dw_v = reshape(Dw,iSys,[]);
Dw_v = Dw_v(:,logicalD);

% interpolate system matrices based on wind speed
A_op_pp = interp1(w_ops,Aw_v,'pchip','pp');
B_op_pp = interp1(w_ops,Bw_v,'pchip','pp');
C_op_pp = interp1(w_ops,Cw_v,'pchip','pp');
D_op_pp = interp1(w_ops,Dw_v,'pchip','pp');

A_op = @(w) ppval(A_op_pp,w);
B_op = @(w) ppval(B_op_pp,w);
C_op = @(w) ppval(C_op_pp,w);
D_op = @(w) ppval(D_op_pp,w);

%% model-dependent elements
% construct interpolated state operating points as a function of wind speed
Xo_pp = pchip(w_ops,x_opsM);
[breaks,coefs,l,k,d] = unmkpp(Xo_pp);

% make the polynomial that describes the derivative
DXo_pp = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
DXo_fun = @(w) ppval(DXo_pp,w);
Xo_fun = @(w) ppval(Xo_pp,w);

% construct interpolated input operating points as a function of wind speed
Uo_pp = spline(w_ops,u_opsM);
Uo_fun = @(w) ppval(Uo_pp,w);

%% scenerio dependent elements
% windcase = 1; %turbulent wind profile
% windcase = 1:2; % step wind profile

for windcase = 1:2

switch windcase
    case 1
        Scenario = load('Turb.mat');        
        idx = contains(Scenario.channame,'Time');
        tt = Scenario.chan(:,idx);
        
        % inputs
        for k = 1:length(indOuts)
            switch indOuts(k)
                case iWindSpeed
                    idx = contains(Scenario.channame,'Wind1VelX');
                    U(:,k) = Scenario.chan(:,idx);
                    iBWindSpeed = k;
                case iGenTorq
                    idx = contains(Scenario.channame,'GenTq');
                    U(:,k) = 1000*Scenario.chan(:,idx);
                    iBGenTorq = k;
                case iBldPitch
                    idx = contains(Scenario.channame,'BldPitch1');
                    U(:,k) = deg2rad(Scenario.chan(:,idx));
                    iBBldPitch = k;
            end
            
            
        end
        
    case 2
        
        Scenario = load('Step.mat');        
        idx = contains(Scenario.channame,'Time');
        tt = Scenario.chan(:,idx);
        
        for k = 1:length(indOuts)
            switch indOuts(k)
                case iWindSpeed
                    idx = contains(Scenario.channame,'Wind1VelX');
                    U(:,k) = Scenario.chan(:,idx);
                    iBWindSpeed = k;
                case iGenTorq
                    idx = contains(Scenario.channame,'GenTq');
                    U(:,k) = 1000*Scenario.chan(:,idx);
                    iBGenTorq = k;
                case iBldPitch
                    idx = contains(Scenario.channame,'BldPitch1');
                    U(:,k) = deg2rad(Scenario.chan(:,idx));
                    iBBldPitch = k;
            end
            
            
        end
        
        
end
%tt = Scenario.Chan.tt;
U_pp = spline(tt,U');
U_fun = @(t) ppval(U_pp,t);

% Min and max values
Wmin = min(U(:,iBWindSpeed));
Wmax = max(U(:,iBWindSpeed));
Gavg = geomean(U(:,iBWindSpeed));

% potentially reduce the number of time points for the inputs
if reducedInputFlag
    N = 50;
    tt = linspace(tt(1),tt(end),N);
    U = U_fun(tt)';
    U_pp = spline(tt,U');
    U_fun = @(t) ppval(U_pp,t);
end

Wind_speed = U(:,iBWindSpeed);

% average wind speed
Wavg = trapz(tt,Wind_speed)/(tt(end)-tt(1));

% construct interpolated wind speed as a function of time
ppW = spline(tt,Wind_speed);
W_fun = griddedInterpolant(tt,Wind_speed,'spline');

[breaks,coefs,l,k,d] = unmkpp(ppW);

% make the polynomial that describes the derivative
ppDW = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
DW_fun = @(t) ppval(ppDW,t);

% construct operating point derivative as a function of time
DXoDt_fun = @(t) (DXo_fun(W_fun(t)))'.*DW_fun(t);

%% optional plots
if debugPlotFlag

    % matrices sparsity
    hf = figure; hf.Color = 'w'; hold on
    subplot(2,2,1)
    spy(Asparsity); title('A');
    subplot(2,2,2)
    spy(Bsparsity); title('B');
    subplot(2,2,3)
    spy(Csparsity); title('C');
    subplot(2,2,4)
    spy(Dsparsity); title('D');

    % time vector
    Tt = linspace(tt(1),tt(end),1e6)';

    % plot interpolated wind speed
    hf = figure; hf.Color = 'w'; hold on
    plot(Tt,W_fun(Tt))
    plot(tt,Wind_speed,'.')

    % plot derivative of interpolated wind speed
    hf = figure; hf.Color = 'w'; hold on
    plot(Tt,DW_fun(Tt))

    % wind speed vector
    Wt = linspace(w_ops(1),w_ops(end),1e4)';

    % plot operating point as a function of wind speed
    hf = figure; hf.Color = 'w'; hold on
    plot(Wt,Xo_fun(Wt))
    plot(w_ops,x_opsM,'.')

    % plot derivative of operating point as a function of wind speed
    hf = figure; hf.Color = 'w'; hold on
    plot(Wt,DXo_fun(Wt))

    % time vector
    tt = linspace(tt(1),tt(end),1e5)';

    % compute operating point derivative as a function of time
    Ws = W_fun(tt);
    DXoDWs = DXo_fun(Ws);
    DwDt = DW_fun(tt);
    DXoDt = DXoDWs'.*DwDt;

    % plot operating point derivative as a function of time
    plot(tt,DXoDt)

end

%--------------------------------------------------------------------------
% time horizon for simulation
TSPAN = [tt(1),tt(end)];
% TSPAN = [tt(1),50]; % shorter simulation

% options
OPTIONS = odeset('reltol',1e-7);

% average wind speed
caseflag = 1;
Y01 =zeros(1,101);Y01(1) = deg2rad(3); Y01(5) = 6/9.5492965964254;Y01(2) = 0; % from ElastoDyn.dat
Y0 = Y01-Xo_fun(Wavg)';
[Tl_avg,Yl_avg] = ode45(@(t,y) odefun(t,y,A_op,B_op,W_fun,U_fun,Uo_fun,DXoDt_fun,Wavg,caseflag,logicalA,logicalB,sizeA,sizeB),...
    TSPAN,Y0,OPTIONS); % run simulation

% wind-varying model
caseflag = 2;
%Y0 =zeros(1,101);Y0(1) = deg2rad(1); Y0(5) = 7.56/9.5492965964254;Y0(2) = 0; % from ElastoDyn.dat
Y0 = Y01 - Xo_fun(W_fun(0))'; % initialize states
[Tl,Yl] = ode45(@(t,y) odefun(t,y,A_op,B_op,W_fun,U_fun,Uo_fun,DXoDt_fun,Wavg,caseflag,logicalA,logicalB,sizeA,sizeB),...
    TSPAN,Y0,OPTIONS); % run simulation

% shift linearized states
Yoff_avg = Yl_avg + Xo_fun(Wavg)';
Yoff = Yl + Xo_fun(W_fun(Tl))';

% shift linearized inputs
Uoff_avg = U_fun(Tl_avg)';
Uoff = U_fun(Tl)';

%% Plot
ComparisonPlot
end
return

%% derivative function
function yp = odefun(t,y,A_op,B_op,W_fun,U_fun,Uo_fun,DXoDt_fun,Wavg,caseflag,logicalA,logicalB,sizeA,sizeB)

switch caseflag
    %----------------------------------------------------------------------
    case 1
        % model wind value
        w = Wavg;

        % state derivative elements
        A = SparseInterp(A_op,logicalA,sizeA,w);
        B = SparseInterp(B_op,logicalB,sizeB,w);
        u = U_fun(t);
        Uo = Uo_fun(w);

        % state derivative function
        yp = A*y + B*(u-Uo);
    %----------------------------------------------------------------------
    case 2
        % model wind value
        w = W_fun(t);

        % state derivative elements
        A = SparseInterp(A_op,logicalA,sizeA,w);
        B = SparseInterp(B_op,logicalB,sizeB,w);
        u = U_fun(t);
        Uo = Uo_fun(w);
        dXoDt = DXoDt_fun(t);

        % state derivative function
        yp = A*y + B*(u-Uo) - dXoDt(:);
    %----------------------------------------------------------------------
end

%disp(num2str(t,'%0.4f'))

end

%% interpolation function with only on non-zero entries
function A = SparseInterp(f,I,s,x)

% initialize
A = zeros(s);

% interpolate
B = f(x);

% assign
A(I) = B;

end