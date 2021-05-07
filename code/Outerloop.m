% Outerloop.m
% Provides the constraints, wind inputs, linearized file to the innerloop
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
function Outerloop

clc; clear;

% pitch contraint
 Constraints = [3,4,5,6]; % list of pitch constraints
Constraints = 6;

% matfile corresponding to the PtfmMass
CaseNames = {'pd_0.2_linear','pd_0.3_linear',...
    'pd_0.4_linear','pd_0.5_linear','pd_0.6_linear',...
    'pd_0.7_linear','pd_0.8_linear','pd_0.9_linear',...
    'pd_1.0_linear','pd_1.1_linear','pd_1.2_linear'};
CaseNames = {'pd_1.0_linear'};

% fraction of nominal mass
Mfactor = [0.2:0.05:0.4,0.5:0.1:1.2];


% Wind Files
WindFiles = {'iea15mw_0.mat','iea15mw_1.mat','iea15mw_2.mat','iea15mw_3.mat','iea15mw_4.mat','072720_183300.mat'};
WindFiles = {'072720_183300.mat'};

% plot
PlotFlag = true;

% DTQP options
opts.dt.nt = 1000;
opts.general.displevel = 2; % 0:none, 1:minimal, 2:verbose

% length
nc = length(Constraints);
nw = length(WindFiles);
nm = length(CaseNames);

% initialize storage
Fflag = false(nm,nc,nw);
StoreVals = cell(nm,nc,nw);

% calculate inner loop response
tic
for m = 1:nm
    for c = 1:nc
        for w = 1:nw
            [T,U,X,Pitch_rate,F] = Innerloop(CaseNames{m},Constraints(c),WindFiles{w},PlotFlag,opts);
            Fflag(m,c,w) = isnan(F);
            StoreVals{m,c,w} = {T,U,X};
        end
    end
end
time = toc;

% LCOE
LCOEflag = 0;

% plot LCOE
LCplot = 1;

if LCOEflag
     [LC,CC,AEP] = LCOE(Fflag,StoreVals,Mfactor,LCplot);
end

end