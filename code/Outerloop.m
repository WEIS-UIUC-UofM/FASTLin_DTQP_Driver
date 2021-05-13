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

%% pitch contraint
 Constraints = [3,4,5,6]; % list of pitch constraints
Constraints = 6;

%% Mass fractions
M = [0.2:0.01:0.3,0.4:0.1:1.2];
M = 1;
% number of masses
nm = length(M);

% initialize
CaseNames = cell(1,nm);
LM = cell(1,nm);
idx = false(1,nm);

for i = 1:nm
    % matfile corresponding to the PtfmMass
    if M(i) == 1
        CaseNames{i} = ['pd_',num2str(M(i),'%.1f'),'_linear.mat'];
    else
        CaseNames{i} = ['pd_',num2str(M(i)),'_linear.mat'];
    end
    
    % check if it exists
    idx(i) = exist(CaseNames{i},'file')==2;
end

% Matfiles that dont exist
ind0 = find(~idx);

% Matfiles that exist
ind1 = find(idx);

% load files tht are available
for i = 1:length(ind1)
    LM{ind1(i)} = load(CaseNames{ind1(i)}) ; 
end


if ~isempty(ind0)
    
    % Find the fractions that dont exist
    for i = 1:length(ind0)
        LinModelFile = CaseNames{ind0(i)};
        ind = regexp(LinModelFile,'\d');
        Mi(i) = str2num(LinModelFile(ind(1):ind(end)));
    end
    
    % Interpolate and find the system matrices for the given fraction
    X = Mass2mat(Mi);
    
    % Assign
    for i = 1:length(ind0)
        LM{ind0(i)} = X{i};
    end
end


%% Wind Files
WindFiles = {'iea15mw_0.mat','iea15mw_1.mat','iea15mw_2.mat','iea15mw_3.mat','iea15mw_4.mat','072720_183300.mat'};
WindFiles = {'072720_183300.mat'};

% plot
PlotFlag = 1;

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
            [T,U,X,Pitch_rate,F] = Innerloop(LM{m},Constraints(c),WindFiles{w},PlotFlag,opts);
            Fflag(m,c,w) = isnan(F);
            StoreVals{m,c,w} = {T,U,X};
        end
    end
end
time = toc;

% LCOE
LCOEflag = 1;

% plot LCOE
LCplot = 1;

if LCOEflag
     [LC,CC,AEP] = LCOE(Fflag,StoreVals,M,LCplot,Constraints);
end

end