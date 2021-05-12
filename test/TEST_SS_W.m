% TEST_SS_W.m
% Test file for interpolating matrices according to wind speed
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
close all; clear; clc;

%pathpdf = mfoldername(mfilename('fullpath'),'S3-FiguresINTMAT');

% load the data
load('pd_1.0_linear.mat')
% test number
testnum = 2;

switch testnum
    %----------------------------------------------------------------------
    case 1
    FitData = 1:20;
    ValidationData = [1];
    %----------------------------------------------------------------------
    case 2
    FitData = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,56];
    ValidationData = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54];
    %----------------------------------------------------------------------
    case 3
    FitData = 1:20;
    ValidationData = [1];
    balflag = true;
    %----------------------------------------------------------------------
    case 4
        
end

%figure

Xop = [SS_Ops_full(:).xop];
Uop = [SS_Ops_full(:).uop];


% important outputs
Io = [1,2,26,12,11,13,14,15,16,41];

% number of fitting and validation models
Nfit = length(FitData);
Nval = length(ValidationData);

% color maps
cfit = [183 28 28]/255;
cint = [13 71 161]/255;
cval = tint(cfit,0.6);

% interpolation method
interpmethod = 'spline';

% transformation matrix (one for all models)
T = eye(size(P{FitData(1)}.A)); % default
Ptb = 1e-3;
Px = zeros(101,1);
Px(1:5) = Ptb;
Pu = Ptb*ones(3,1);

% create data matrices
for k = 1:Nfit    
    A(k,:,:) = (T*P{FitData(k)}.A)/T;
    Ao(k,:,:) = P{FitData(k)}.A;
    B(k,:,:) = T*P{FitData(k)}.B;
    Bo(k,:,:) = P{FitData(k)}.B;
    C(k,:,:) = P{FitData(k)}.C/T;
    D(k,:,:) = P{FitData(k)}.D;
    X(k,:) = Xop(:,k);
    U(k,:) = Uop(:,k);
    Dx(k,:) = squeeze(Ao(k,:,:))*(X(k,:)'+Px) + squeeze(Bo(k,:,:))*(U(k,:)'+Pu);
end

% current wind speed values
W = WindSpeed(FitData);

% 
Winterp = linspace(WindSpeed(1),WindSpeed(end),1000);

%% A matrix
figure('Color','w'); hold on

% interpolate matrix
Ainterp = interp1(W,A,Winterp,interpmethod);

% permute dimensions (so squeeze isn't needed)
Ainterp = permute(Ainterp,[2,3,1]);

% plot interpolated data points
ha1 = plotA(Ainterp,Winterp,cint,6);

% plot fitting data points
for k = 1:Nfit    
    Afit(:,:,k) = (T*P{FitData(k)}.A)/T;
   ha2 = plotA(Afit(:,:,k),WindSpeed(FitData(k)),cfit,16);
end

% plot validation points
for k = 1:Nval
    Aval(:,:,k) = (T*P{ValidationData(k)}.A)/T;
   ha3 = plotA(Aval(:,:,k),WindSpeed(ValidationData(k)),cval,16);
end

xlabel('Real','FontSize',14);xlim([-0.1653,-0.1638])
ylabel('Imag','FontSize',14); ylim([0.74,0.7415])
legend([ha1 ha2 ha3],{'Interpolated Points','Fitting Data','Validation Data'},'FontSize',12,'Orientation','horizontal','location','northoutside','Interpreter','latex');

%title('Pole Zero Map')

%% B matrix
figure
ha = gca; ha.LineWidth = 1; ha.FontSize = 14;
hf = gcf; hf.Color = 'w';
hold on;

% interpolate matrix
Binterp = interp1(W,B,Winterp,interpmethod);

% permute dimensions (so squeeze isn't needed)
Binterp = permute(Binterp,[2,3,1]);

% normalize
Bmax = max(abs(Binterp),[],3) + 1e-4;

% plot interpolated data points
%plotBCD(Binterp./Bmax,Winterp,'b',6,'-')
Bx = Binterp./Bmax;

hb1 = plot(Winterp,reshape(Bx(:),[],length(Winterp))','color',cint,'markersize',6);

% plot fitting data points
for k = 1:Nfit
    Bfit(:,:,k) = T*P{FitData(k)}.B;
end
%plotBCD(Bfit./Bmax,WindSpeed(FitData),'k',16,'.')
Bx = Bfit./Bmax;
hb2 = plot(WindSpeed(FitData),reshape(Bx(:),[],length(WindSpeed(FitData)))','.','color',cfit,'markersize',16);

% plot validation points
for k = 1:Nval
    Bval(:,:,k) = T*P{ValidationData(k)}.B;
end
%plotBCD(Bval./Bmax,WindSpeed(ValidationData),'r',16,'.')
Bx = Bval./Bmax;

hb3 = plot(WindSpeed(ValidationData),reshape(Bx(:),[],length(WindSpeed(ValidationData)))','.','color',cval,'markersize',16);

% y label for plot
ylabel('Normalized B matrix entries','FontSize',14)
xlabel('Wind Speed','FontSize',14);xlim([3,25])
legend([hb1(1),hb2(2),hb3(3)],{'Interpolated Points','Fitting Points','Validation Points'},'FontSize',12,'Orientation','horizontal','location','northoutside','Interpreter','latex');

%% operating points

figure
ha = gca; ha.LineWidth = 1; ha.FontSize = 14;
hf = gcf; hf.Color = 'w';
hold on;
% interpolate

Xinterp = interp1(W,X,Winterp,interpmethod);
Dxinterp = interp1(W,Dx,Winterp,interpmethod);

% normalize
Dxmax = max(abs(Dxinterp'),[],2) + 1e-4;

% wind speed corresponding to data indices
Valind = WindSpeed(ValidationData);
Fitind = WindSpeed(FitData);

Dact = zeros(Nfit,1);
Dval = zeros(Nval,1);

for k = 1:Nfit
   Xfit(k,:) = Xop(:,k);
   Ufit(k,:) = Uop(:,k);
   Dxfit(k,:) = squeeze(Afit(:,:,k))*(Xfit(k,:)'+Px) + squeeze(Bfit(:,:,k))*(Ufit(k,:)'+Pu);
   Dact(k) = norm(Dx(k,:)-Dxfit(k,:),inf);
end



for k = 1:Nval
   Xval(k,:) = Xop(:,k); 
   Uval(k,:) = Uop(:,k);
   Dxval(k,:) = squeeze(Aval(:,:,k))*(Xval(k,:)'+Px) + squeeze(Bval(:,:,k))*(Uval(k,:)'+Pu);
   Dval(k) = norm(Dx(k,:)-Dxval(k,:),inf);
end

% plot
hx1 = plot(Winterp,Xinterp,'Color',cint,'markersize',6);
hx2 = plot(WindSpeed(FitData),Xfit,'.','color',cfit,'markersize',16);
hx3 = plot(WindSpeed(ValidationData),Xval,'.','color',cval,'markersize',16);
ylabel('State Operating points values','FontSize',14)
xlabel('Wind Speed','FontSize',14);xlim([3,25])
legend([hx1(1) hx2(2) hx3(3)],{'Interpolated Points','Fitting Data','Validation Data'},...
    'FontSize',12,'Orientation','horizontal','location','northoutside','Interpreter','latex');


convertflag = 0;

Dxfit = Dxfit./Dxmax';
Dxval = Dxval./Dxmax';
% Dxinterp = Dxinterp./Dxmax';

Wcombined = [Valind,Fitind];
Dcombined = [Dval',Dact'];
[Wcombined,Isort] = sort(Wcombined);
Dcombined = Dcombined(Isort);

figure
ha = gca; ha.LineWidth = 1; ha.FontSize = 14;
hf = gcf; hf.Color = 'w';
hold on;
hx1 = plot(WindSpeed(FitData),Dact,'.','color',cfit,'markersize',16);
hx2 = plot(WindSpeed(ValidationData),Dval,'.','color',cval,'markersize',16);
hx3 =  plot(Wcombined,Dcombined,'-','Color',cint,'linewidth',1.5);

ylabel('$H_{\infty}$ of ($\dot{x}_{act}-\dot{x}_{int}$)','FontSize',14,'Interpreter','Latex')
xlabel('Wind Speed','FontSize',14);xlim([3,25])
hx3 = legend([hx1 hx2],'Training Points','Validation Points','FontSize',14,'Orientation','vertical','location','best');

%%

if convertflag
    
    ifx = findobj('type','figure');
    
    n = length(ifx);
    
    Savenames = {'EigA','IntpB','Xoper'};
    
    for i = 1:n
        
        figure(i)
        filename = [pathpdf,Savenames{i}];
        str = ['export_fig ''',filename,''' -pdf'];
        eval(str)
        
    end
    
end

%%
% plot matrix values
function plotBCD(B,W,c,varargin)

if ~isempty(varargin)
	markersize = varargin{1};
    linestyle = varargin{2};
else
	markersize = 16;
    linestyle = '.';
end

% plot matrix
plot(W,reshape(B(:),[],length(W))',linestyle,'markersize',markersize,'color',c)
%legend('Interpolated Points','Fitting Data','Validation Data','FontSize',14,'Orientation','horizontal','location','northoutside','Interpreter','latex');
end

function PlotX(X,W,c,varargin)
if ~isempty(varargin)
	markersize = varargin{1};
    linestyle = varargin{2};
else
	markersize = 16;
    linestyle = '.';
end

plot(W,X,linestyle,'markersize',markersize,'color',c)

end

% plot square matrix eigenvalues
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