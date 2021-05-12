%--------------------------------------------------------------------------
% ComparisonPlot.m
% Plot file for TEST_Time_Dom.m
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
% load data

mname = mfilename('fullpath');
[mpath,mname] = fileparts(mname);
savepath = fullfile(mpath,'SaveData',mname);

Convertflag = false;

%% Plot
%% Wind speed
figure(1)
ha = gca; ha.LineWidth = 1; ha.FontSize = 14;
hf = gcf; hf.Color = 'w';
%subplot(2,1,i)
idx = contains(Scenario.channame,'Wind1VelX');
plot(Scenario.chan(:,1),Scenario.chan(:,idx),'linewidth',1,'markersize',1);hold on;
xlabel('Time [s]','FontSize',14,'Interpreter','latex');ylabel('Wind Speed [m/s]','FontSize',14,'Interpreter','latex');
legend('Study 1','Study 2','FontSize',14,'Orientation','horizontal','location','northoutside');


%% PtfmPitch

figure
ha = gca; ha.LineWidth = 1; ha.FontSize = 14;
hf = gcf; hf.Color = 'w';
hold on;
% Linearized
plot(Tl_avg,rad2deg(Yoff_avg(:,iPltPitch)),'-','linewidth',2,'markersize',1);
plot(Tl,rad2deg(Yoff(:,iPltPitch)),'-','linewidth',2,'markersize',1);
% Nonlinear
idx = contains(Scenario.channame,'PtfmPitch');
plot(Scenario.chan(:,1),Scenario.chan(:,idx),'-','linewidth',2,'markersize',1);
xlabel('Time [s]','FontSize',14,'Interpreter','latex');ylabel('Platform Tilt [deg]','FontSize',14,'Interpreter','latex');

legend('$w_{avg}$','$w(t)$','OpenFAST','FontSize',14,'Orientation','horizontal','location','northoutside','Interpreter','latex');

%% GenSpeed
 figure
 ha = gca; ha.LineWidth = 1; ha.FontSize = 14;
 hf = gcf; hf.Color = 'w';
 hold on
% Linearized
plot(Tl_avg,Yoff_avg(:,iGenSpeed),'-','linewidth',2,'markersize',1);
plot(Tl,Yoff(:,iGenSpeed),'-','linewidth',2,'markersize',1);
% Nonlinear
idx = contains(Scenario.channame,'GenSpeed');
plot(Scenario.chan(:,1),(Scenario.chan(:,idx))/9.5492965964254,'-','linewidth',2,'markersize',1);
xlabel('Time [s]','FontSize',14,'Interpreter','latex');ylabel('GenSpeed [rad/s]','FontSize',14,'Interpreter','latex');
legend('$w_{avg}$','$w(t)$','OpenFAST','FontSize',14,'Orientation','horizontal','location','northoutside','Interpreter','latex');


%% Convert

if Convertflag
    %addpath(genpath('A:\832592734\IDETC-code\export_fig'));
    
    pathpdf = mfoldername(mfilename('fullpath'),'S3-Figures');
    
    ifx = findobj('type','figure');
    
    n = length(ifx);
    
    Savenames = {'WindInputs','S1PtfmPitch','S1GenSpeed','S2PtfmPitch','S2GenSpeed'};
    
    for i = 1:n
        
        figure(i)
        filename = [pathpdf,Savenames{i}];
        str = ['export_fig ''',filename,''' -pdf'];
        eval(str)
        
    end
end