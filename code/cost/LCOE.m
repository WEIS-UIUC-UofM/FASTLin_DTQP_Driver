% LCOE.m
% Calculates and plots the LCOE values for the different mass and
% constraint cases
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
function [LC,CC,AEP] = LCOE(Fflag,StoreVals,Mfactor,LCplot,c)
% determine size

s = size(Fflag,[1 2 3]);
nm = s(1);nc = s(2);nw = s(3);

LD = cell(1,nc);

for i = 1:nc
    LD{i} =  ['$\Theta_p$ $\leq$',num2str(c(i)),' deg'];
end

% Calculate system cost
[Mturbine,Cturbine] = SystemCost(Mfactor);

% Tunable parameters

% number of years
n = 30;
k = 1:n;

% discount rate
r = 0.07;
NormFactor = 1./(1+r).^k;

% initialize
LC = zeros(nm,nc);
CC = zeros(nm,nc);
AEP = zeros(nm,nc);
CapFactor = zeros(nm,nc);

% weights for the DLCs
Wts = [0.0631,0.2525,0.3156,0.1479,0.1262,0.0947];


for m =  1:nm
    for  c = 1:nc
        for w = 1:nw
            if ~Fflag(m,c,w) == 1 % check for infeasible values
                % extract data
                ILdata = StoreVals{m,c,w};
              
                % states and controls
 
                U = ILdata{2};
                X = ILdata{3};

                % Combine the wind speed and genpwr
                Ws(:,w) = U(:,1);
                GenPwr(w) = mean(U(:,2).*X(:,5)*0.9655/10^6);
            else
                GenPwr(w) = 0;
                
            end
        end
        try
            % downtime
            GenPwr = 0.85*GenPwr;

            % AEP
            E = sum(Wts.*GenPwr);

            % Capacity factor
            CapFactor(m,c) = E/15;

            % production over a year
            E = E*(24*365);
        catch
            E = 0;
        end

         % En
        Den = NormFactor.*E;
        Den = sum(Den);
        
        % Cn
        Op = 0.02*Den;
        Num = Cturbine(m)/(1+r) + Op;
        
        % Cn/En
        LC(m,c) = Num/Den;
        CC(m,c) = Num;
        AEP(m,c) = Den;
    end
    Ws = [];
    GenPwr = [];
    ILdata = [];
end

% plot flag

if LCplot
    
    
    figure
    ha = gca; ha.LineWidth = 1; ha.FontSize = 14; %ha.TickLabelInterpreter='latex';
    hf = gcf; hf.Color = 'w';
    
    plot(Mfactor*100,LC,'.-','LineWidth',2,'MarkerSize',15)
    xlabel('Percentage of Nominal Platform Mass [\%]','FontSize',14,'Interpreter','latex');
    ylabel('LCOE [\$/MWh]','FontSize',14,'Interpreter','latex');
    Ld = legend(LD,...
        'FontSize',14,'location','northoutside','Interpreter','latex','Orientation','horizontal');
    Ld.NumColumns = 3;
end
end