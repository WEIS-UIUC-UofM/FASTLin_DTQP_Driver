% TEST_Freq_Dom.m
% Frequency domain comparison between interpolated and actual models
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
clc; clear; close all;
LinearModels = load('pd_1.0_linear.mat');

original_color = [183 28 28]/255;
interpolated_color = [13 71 161]/255;
validation_color = tint(original_color,0.6);
bcolor = [0 0 0]; % black

%% Extract

nl = length(LinearModels.P);

FitData = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,56];
ValidationData = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54];

%pathpdf = mfoldername(mfilename('fullpath'),'S3-FiguresFD');

% length of data matrices
Nfit = length(FitData);
Nval = length(ValidationData);

% go through each linear model
for iSys = 1:Nfit
    
    k = FitData(iSys);
    
    % extract
    sys = LinearModels.P{k};
    
    % number of states
    nStates(iSys) = size(LinearModels.P{k}.A,1);
    
    % wind speed
    w_ops(iSys,1) = LinearModels.WindSpeed(k);
    
    % matrices
    Am(:,:,iSys) = LinearModels.P{k}.A;
    Bm(:,:,iSys) = LinearModels.P{k}.B;
    Cm(:,:,iSys) = LinearModels.P{k}.C;
    Dm(:,:,iSys) = LinearModels.P{k}.D;
    
end
% permute matrices
Am = permute(Am,[3,1,2]);
Bm = permute(Bm,[3,1,2]);
Cm = permute(Cm,[3,1,2]);
Dm = permute(Dm,[3,1,2]);

% interpolate system matrices based on wind speed


% wind speed corresponding to data indices
Valind = LinearModels.WindSpeed(ValidationData);
Fitind = LinearModels.WindSpeed(FitData);

% output
OUT_ind = find(contains(LinearModels.P{1}.OutputName,'ED GenSpeed'));% 'ED PtfmPitch' 'ED GenSpeed'
out_name = LinearModels.P{1}.OutputName(OUT_ind);

% input
IN_ind = find(contains(LinearModels.P{1}.InputName,'ED BldPitchCommand'));% 'ED PtfmPitch' 'ED GenSpeed'
in_name = LinearModels.P{1}.InputName(IN_ind);

% initalize
Nrm = zeros(1,Nval);
act_cell = cell(Nval,1);
int_cell = cell(Nval,1);

sys_intp = cell(Nval,1);
sys_actp = cell(Nval,1);

%interp_method = {'pchip'};

% interpolate
A_op = @(w) interp1(w_ops,Am,w,'pchip');
B_op = @(w) interp1(w_ops,Bm,w,'pchip');
C_op = @(w) interp1(w_ops,Cm,w,'pchip');
D_op = @(w) interp1(w_ops,Dm,w,'pchip');

for i = 1:Nval
    
    % get index
    idx = ValidationData(i);
    
    % find actual matrix vlue
    A_act = LinearModels.P{idx}.A;
    B_act = LinearModels.P{idx}.B;
    C_act = LinearModels.P{idx}.C;
    D_act = LinearModels.P{idx}.D;
    
    % build ss
    sys_act = ss(A_act,B_act(:,IN_ind),C_act(OUT_ind,:),D_act(OUT_ind,IN_ind));
    sys_actp{i} = sys_act;
    % set up transfer function
%     [b_actp,a_actp] = ss2tf(sys_act.A,sys_act.B,sys_act.C,sys_act.D,3);
    
    % store value
%     sys_actp{i} = tf(b_actp,a_actp);
    
    % interpolated matrix
    ws = Valind(i);
    
    % interpolated value
    A_int = squeeze(A_op(ws));
    B_int = squeeze(B_op(ws));
    C_int = squeeze(C_op(ws));
    D_int = squeeze(D_op(ws));
    
    % build ss
    sys_int = ss(A_int,B_int(:,IN_ind),C_int(OUT_ind,:),D_int(OUT_ind,IN_ind));
    sys_intp{i} = sys_int;
    % set up transfer function
%     [b_intp,a_intp] = ss2tf(sys_int.A,sys_int.B,sys_int.C,sys_int.D,3);
    % store value
%     sys_intp{i} = tf(b_intp,a_intp);
    
    Nrm(i) = norm((sys_actp{i}-sys_intp{i}),inf);
    % Hinf norm of the systems
    Nint(i) = norm(sys_act-sys_int,inf);
    
    % store value
    act_cell{i} = sys_act;
    int_cell{i} = sys_int;
end

Nact = zeros(1,Nfit);

% fitting data
for i = 1:Nfit
    
    idx = FitData(i);
    
    % original system
    A_act = LinearModels.P{idx}.A;
    B_act = LinearModels.P{idx}.B;
    C_act = LinearModels.P{idx}.C;
    D_act = LinearModels.P{idx}.D;
    
    
    % build ss
    sys_act = ss(A_act,B_act,C_act(OUT_ind,:),D_act(OUT_ind,:));
    
    % interpolated matrix
    ws = Fitind(i);
    
    % interpolate the system
    A_int = squeeze(A_op(ws));
    B_int = squeeze(B_op(ws));
    C_int = squeeze(C_op(ws));
    D_int = squeeze(D_op(ws));
    
    % build ss
    sys_int = ss(A_int,B_int,C_int(OUT_ind,:),D_int(OUT_ind,:));
    
    % find norm
    Nact(i) = norm(sys_act-sys_int,inf);
    
    
end

%%
figure
hold on;
ha = gca; ha.LineWidth = 1; ha.FontSize = 14;
hf = gcf; hf.Color = 'w';

% bode plot
set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter\

% combine data
Wcombined = [Valind,Fitind];
Hinfcombined = [Nint,Nact];
[Wcombined,Isort] = sort(Wcombined);
Hinfcombined = Hinfcombined(Isort);

% plot interpolated and actual values
h1 = plot(Wcombined,Hinfcombined,'-','Color',validation_color,'linewidth',1.5);
h3 = plot(Valind,Nint,'.','Color',validation_color,'markersize',20);
h2 = plot(Fitind,Nact,'.','Color',original_color,'markersize',20);

hl = legend([h2 h3],'Training Points','Validation Points','FontSize',14,'Orientation','vertical','location','best');
hl.EdgeColor = bcolor;

ylabel('$H_{\infty}$ norm Error');
xlabel('Wind Speed [m/s]');

ha = gca; ha.LineWidth = 1; ha.FontSize = 16;
ha.GridColor = bcolor;
ha.XColor = bcolor;
ha.YColor = bcolor;
ha.MinorGridColor = bcolor;
% ha.XScale = 'log';
% ha.YScale = 'linear';
xlim([min(Wcombined) max(Wcombined)])


%%
figure
hf = gcf; hf.Color = 'w'; hold on

% find worst point
[M,I] = max(Nint);

disp(['max error at w=',num2str(Valind(I)),'m/s'])
% extract system
actp = sys_actp{I};
intp = sys_intp{I};
n_max = norm(actp-intp,inf);

% bode plot
set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter

W = logspace(-3,2,1e4);
[MAG1,PHASE1] = bode(actp,W);
[MAG2,PHASE2] = bode(intp,W);

plot(W,mag2db(squeeze(MAG1)),'linewidth',2,'Color',original_color)
plot(W,mag2db(squeeze(MAG2)),'linewidth',2,'Color',interpolated_color)

hl = legend('Original System','Interpolated System','FontSize',14,'Orientation','vertical','location','best');
hl.EdgeColor = bcolor;

ha = gca; ha.LineWidth = 1; ha.FontSize = 16;
ha.GridColor = bcolor;
ha.XColor = bcolor;
ha.YColor = bcolor;
ha.MinorGridColor = bcolor;
ha.XScale = 'log';
ha.YScale = 'linear';
xlim([min(W) max(W)])
xlabel('$\omega$ [rad/s]')
ylabel('Magnitude [dB]')

%%
figure
hf = gcf; hf.Color = 'w'; hold on

% find best point
[M,I] = min(Nint);
disp(['min error at w =',num2str(Valind(I)),'m/s'])

sys_act = act_cell{I};
sys_int = int_cell{I};

% extract system
actp = sys_actp{I};
intp = sys_intp{I};
n_min = norm(actp-intp,inf);
% bode plot
set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter

W = logspace(-3,2,1e4);
[MAG1,PHASE1] = bode(actp,W);
[MAG2,PHASE2] = bode(intp,W);

plot(W,mag2db(squeeze(MAG1)),'linewidth',2,'Color',original_color)
plot(W,mag2db(squeeze(MAG2)),'linewidth',2,'Color',interpolated_color)

hl = legend('Original System','Interpolated System','FontSize',14,'Orientation','vertical','location','best');
hl.EdgeColor = bcolor;

ha = gca; ha.LineWidth = 1; ha.FontSize = 16;
ha.GridColor = bcolor;
ha.XColor = bcolor;
ha.YColor = bcolor;
ha.MinorGridColor = bcolor;
ha.XScale = 'log';
ha.YScale = 'linear';
xlim([min(W) max(W)])
xlabel('$\omega$ [rad/s]')
ylabel('Magnitude [dB]')

%% export options
exportflag = 0;

if exportflag
    Savenames{1} = 'TF_Hinf';
    Savenames{2} = 'TF_maxerror';
    Savenames{3} = 'TF_minerror';
    
    figure(1)
    filename = [pathpdf,Savenames{1}];
    str = ['export_fig ''',filename,''' -pdf'];
    eval(str)
    
    figure(2)
    filename = [pathpdf,Savenames{2}];
    str = ['export_fig ''',filename,''' -pdf'];
    eval(str)
    
    figure(3)
    filename = [pathpdf,Savenames{3}];
    str = ['export_fig ''',filename,''' -pdf'];
    eval(str)
    
end