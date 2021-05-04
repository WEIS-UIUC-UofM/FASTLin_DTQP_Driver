clc; clear; close all;

% Mass cases
CaseNames = {'pd_0.2_linear','pd_0.3_linear','pd_0.4_linear','pd_0.5_linear','pd_0.6_linear',...
    'pd_0.7_linear','pd_0.8_linear','pd_0.9_linear',...
    'pd_1.0_linear','pd_1.1_linear','pd_1.2_linear'};

% Validation and fitting indices
Fitind = 1:2:11;Nfit = length(Fitind);
Valind = 2:2:10;Nval = length(Valind);
interpmethod = 'spline';

% mass values
Mfit = 0.2:0.2:1.2;
Mval = 0.3:0.2:1.1 ;

% Initialize
FitCases = cell(1,Nfit);
ValCases = cell(1,Nval);

% store values
for i = 1:Nfit
FitCases{i} = CaseNames{Fitind(i)};
end

for i = 1:Nval
ValCases{i} = CaseNames{Valind(i)};
end

% initialize
Afit = zeros(Nfit,46,101,101);
Bfit = zeros(Nfit,46,101,3);
Xfit = zeros(Nfit,46,101,1);
Ufit = zeros(Nfit,46,3,1);

for i = 1:Nfit
    LinModelFile = strcat(FitCases{i},'.mat');
    LinearModels = load(LinModelFile);
    
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
        y_opsM(:,iSys) = LinearModels.SS_Ops(iSys).yop;
        
        % matrices
        Am(:,:,iSys) = LinearModels.P{iSys}.A;
        Bm(:,:,iSys) = LinearModels.P{iSys}.B;
        Cm(:,:,iSys) = LinearModels.P{iSys}.C;
        Dm(:,:,iSys) = LinearModels.P{iSys}.D;
        
    end
    
    % permute matrices
    Afit(i,:,:,:) = permute(Am,[3,1,2]);
    Bfit(i,:,:,:) = permute(Bm,[3,1,2]);
    Xfit(i,:,:) = x_opsM';
    Ufit(i,:,:) = u_opsM';
    

    
end

    Xmax = max(x_opsM,[],2) + 1e-4;
    Umax = max(u_opsM,[],2) + 1e-4;
% size of A
sizeA = size(Afit,3,4);nx = sizeA(1);Xind = 1:nx;

% size of B
sizeB = size(Bfit,3,4);nu = sizeB(2);Uind = 1:nu;

% number of wind speed
sizeW = length(w_ops);

% index for specific wind speed
Navg = find(w_ops==11);

% interpolate
Aval = interp1(Mfit,Afit,Mval,'spline');
Bval = interp1(Mfit,Bfit,Mval,'spline');
Xval = interp1(Mfit,Xfit,Mval,'spline');
Uval = interp1(Mfit,Ufit,Mval,'spline');

%% Validation
Valflag = 1;

if Valflag
    % current wind speed values
    M = Mfit;
    
    %
    Minterp = linspace(Mfit(1),Mfit(end),1000);
    
    %% A matrix
    figure('Color','w'); hold on
    
    T = eye(101);
    % color maps
    cfit = spring(Nfit);
    cval = winter(Nval);
    
    % interpolate matrix
    Ainterp = interp1(M,Afit,Minterp,interpmethod);
    
    Ainterp = squeeze(Ainterp(:,Navg,:,:));
    % permute dimensions (so squeeze isn't needed)
    Ainterp = permute(Ainterp,[2,3,1]);
    
    % plot interpolated data points
    ha1 = plotA(Ainterp,Minterp,'k',6);
    
    % plot fitting data points
    for k = 1:Nfit
        ha2 = plotA(squeeze(Afit(k,Navg,:,:)),Mfit(k),cfit(k,:),16);
    end
    
    % plot validation points
    for k = 1:Nval
        ha3 = plotA(squeeze(Aval(k,Navg,:,:)),Mval(k),cval(k,:),16);
    end
    
    xlabel('Real','FontSize',14);xlim([-0.09 -0.045])
    ylabel('Imag','FontSize',14); ylim([-0.6 0.6])
    legend([ha1 ha2 ha3],{'Interpolated Points','Fitting Data','Validation Data'},'FontSize',12,'Orientation','horizontal','location','northoutside','Interpreter','latex');
    
    % plot
    
    figure
    ha = gca; ha.LineWidth = 1; ha.FontSize = 14;
    hf = gcf; hf.Color = 'w';
    hold on;
    
    Xinterp = interp1(M,Xfit,Minterp,interpmethod);
    
    Xinterp = squeeze(Xinterp(:,Navg,:));
    
    xmax = max(abs(Xinterp'),[],2) + 1e-4;
    
    hx1 = plot(Minterp,Xinterp,'b-','markersize',6);
    hx2 = plot(Mfit,squeeze(Xfit(:,Navg,:)),'k.','markersize',16);
    hx3 = plot(Mval,squeeze(Xval(:,Navg,:)),'r.','markersize',16);
    
    ylabel('State Operating points values','FontSize',14)
    xlabel('Mass Fraction','FontSize',14);
    legend([hx1(1) hx2(2) hx3(3)],{'Interpolated Points','Fitting Data','Validation Data'},...
        'FontSize',12,'Orientation','horizontal','location','northoutside','Interpreter','latex');

end

%% Convert 
Cflag = 0;

if Cflag
    P = cell(1,sizeW);
    C = zeros(41,101);
    D = zeros(41,3);
    WindSpeed = w_ops;
    
    for i = 1: length(Mval)
        for j = 1:sizeW
            
            A = squeeze(Aval(i,j,:,:));
            B = squeeze(Bval(i,j,:,:));
            
            P{j} = ss(A,B,C,D);
            SS_Ops(j).xop = squeeze(Xval(i,j,:));
            SS_Ops(j).uop = squeeze(Uval(i,j,:));
        end
        savename = strcat('pd_',num2str(Mval(i)),'_linear.mat');
        save(savename,'P','SS_Ops','WindSpeed')
    end
end

%%

function h = plotA(A,W,c,varargin)

if ~isempty(varargin)
	markersize = varargin{1};
else
	markersize = 16;
end

for k = 1:length(W)
    % eigenvalues of A matrix
    Aeig = eig(A(:,:,k));

    % plot
   h = plot(real(Aeig),imag(Aeig),'.','markersize',markersize,'color',c);

end

end

