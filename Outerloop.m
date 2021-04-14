function Outerloop

% pitch contraint
Constraints = 6;
    
% matfile corresponding to the mass
%CaseNames = {'pd_0.2_linear','pd_0.3_linear','pd_0.4_linear','pd_0.5_linear','pd_0.6_linear',...
%    'pd_0.7_linear','pd_0.8_linear','pd_0.9_linear',...
%    'pd_1.0_linear','pd_1.1_linear','pd_1.2_linear'};
CaseNames = {'pd_0.2_linear'};

% Wind Files
%WindFiles = {'\iea15mw_0.mat','\iea15mw_1.mat','\iea15mw_2.mat','\iea15mw_3.mat','\iea15mw_4.mat','072720_183300.mat'};
WindFiles = {'\iea15mw_1.mat'};

% get size
nm = length(CaseNames);
nw = length(WindFiles);
nc = length(Constraints);

% initialize storage
StoreValues = cell(nm,nw,nc);


for m = 1:nm
    for w = 1:nw
        for c = 1:nc
        
        % Plant file    
        casefile = CaseNames{m};
        
        % Wind file
        WindFile = WindFiles{w};
        
        % Constraints
        PitchConstraint = Constraints(c);
        
        % calculate inner loop response
        [T1,U,X,Pitch_rate,F] = Innerloop(casefile,PitchConstraint,WindFile);
        
        % store
        StoreValues{m,w,c} = {T1,U,X,F};
        
        end
    end        
end 

end