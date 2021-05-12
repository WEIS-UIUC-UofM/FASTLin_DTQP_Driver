% Mass2mat.m
% Interpolates and finds system matrices for intermediate mass fractions
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
function LM = Mass2mat(M)

LM = cell(1,length(M));

% Mass cases
CaseNames = {'pd_0.2_linear','pd_0.3_linear','pd_0.4_linear','pd_0.5_linear','pd_0.6_linear',...
    'pd_0.7_linear','pd_0.8_linear','pd_0.9_linear',...
    'pd_1.0_linear','pd_1.1_linear','pd_1.2_linear'};

% Validation and fitting indices
Fitind = 1:1:11;Nfit = length(Fitind);
interpmethod = 'spline';

% mass values
Mfit = [0.2:0.1:1.2];

% Initialize
FitCases = cell(1,Nfit);


% store values
for i = 1:Nfit
    FitCases{i} = CaseNames{Fitind(i)};
end

% Load files
for i = 1:Nfit
    LinModelFile = strcat(FitCases{i},'.mat');
    LinearModel = load(LinModelFile);
    
    % Set length
    nl = length(LinearModel.P);
    
    % go through each linear model
    for iSys = 1:nl
        
        % extract
        sys = LinearModel.P{iSys};
        
        % number of states
        nStates(iSys) = size(LinearModel.P{iSys}.A,1);
        
        % wind speed
        w_ops(iSys,1) = LinearModel.WindSpeed(iSys);
        
        % operating points
        x_opsM(:,iSys) = LinearModel.SS_Ops(iSys).xop;
        u_opsM(:,iSys) = LinearModel.SS_Ops(iSys).uop;
        
        % matrices
        Am(:,:,iSys) = LinearModel.P{iSys}.A;
        Bm(:,:,iSys) = LinearModel.P{iSys}.B;
        
    end
    
    % permute matrices
    Afit(i,:,:,:) = permute(Am,[3,1,2]);
    Bfit(i,:,:,:) = permute(Bm,[3,1,2]);
    Xfit(i,:,:) = x_opsM';
    Ufit(i,:,:) = u_opsM';
   
end

% number of wind speed
sizeW = length(w_ops);

%% interpolate
Aval = interp1(Mfit,Afit,M,interpmethod);
Bval = interp1(Mfit,Bfit,M,interpmethod);
Xval = interp1(Mfit,Xfit,M,interpmethod);
Uval = interp1(Mfit,Ufit,M,interpmethod);

%% Convert


P = cell(1,sizeW);
C = zeros(41,101);
D = zeros(41,3);
WindSpeed = w_ops;

for i = 1: length(M)
    for j = 1:sizeW
        
        A = squeeze(Aval(i,j,:,:));
        B = squeeze(Bval(i,j,:,:));
        
        P{j} = ss(A,B,C,D);
        SS_Ops(j).xop = squeeze(Xval(i,j,:));
        SS_Ops(j).uop = squeeze(Uval(i,j,:));
    end
    LinearModels.P = P;
    LinearModels.WindSpeed = WindSpeed';
    LinearModels.SS_Ops = SS_Ops;
    LM{i} = LinearModels;
end

end