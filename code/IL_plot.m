% IL_plot.m
% Innerloop plot code
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
close all;

ha = gca; ha.LineWidth = 1; ha.FontSize = 14;
hf = gcf; hf.Color = 'w';
set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter\

% plot

subplot(3,2,1);
plot(T1,U(:,1),'LineWidth',1);
title('Horizontal Wind Speed [m/s]','FontSize',8);

subplot(3,2,5);
plot(T1,U(:,2)/1000,'LineWidth',1);
title('Generator Torque [KN-m]','FontSize',8);
%ylim([1.8e+04 2e+04]);

subplot(3,2,3);
plot(T1,(U(:,3)),'LineWidth',1);
title('Blade Pitching Control [rad]','FontSize',8);

subplot(3,2,4);
plot(T1,0.9655*X(:,5).*U(:,2)/1000,'LineWidth',1);
title('Gen speed * Gen Torque [kW]','FontSize',8);

subplot(3,2,6);
plot(T1,X(:,5),'LineWidth',1);
title('Generator Speed [rad/s]','FontSize',8);
%yticks([0.75 0.775 0.8]);

subplot(3,2,2);
tmp2 = plot(T1,rad2deg(X(:,1)),'LineWidth',1);
title('Platform Tilt [deg]','FontSize',8);
