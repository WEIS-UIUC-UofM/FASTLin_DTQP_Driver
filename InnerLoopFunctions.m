% InnerLoopFunctions.m
% Contains a list of commands to calculate various quantities and perform
% specific functions
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
switch flag
    
    case 'Disc2Cont'
        % This script interpolates the discrete values of the state and control
        % operating points into continuous functions dependednt on wind speed
        % Specific functions of GeneratorSpeed, GeneratorTorque and BladePitch are
        % required in the calculation of the objective term
        
        %% Xo_fun
        % interpolate discrete state operating points into a continuous function on
        % 'w'
        
        % State operating points
        Xo_pp = pchip(w_ops,x_opsM);
        Xo_fun = @(w) ppval(Xo_pp,w);
        
        % Generator Speed
        GS_pp = pchip(w_ops,x_opsM(5,:));
        GS_fun = @(w)ppval(GS_pp,w);
        
        % Generator Speed*-1
        GSn_pp = pchip(w_ops,-x_opsM(5,:));
        GSn_fun = @(w)ppval(GSn_pp,w);
        
        %% Uo_fun
        % interpolate discrete control operating points into a continuous function
        % on 'w'
        
        % Control operating points
        Uo_pp = pchip(w_ops,u_opsM);
        Uo_fun = @(w) ppval(Uo_pp,w);
        
        % Generator Torque
        GenTq_pp = pchip(w_ops,u_opsM(2,:));
        GT_fun = @(w) ppval(GenTq_pp,w);
        
        % negative Generator Torque
        GenTqn_pp = pchip(w_ops,-u_opsM(2,:));
        GTn_fun = @(w) ppval(GenTqn_pp,w) ;
        
        % Blade Pitch
        BPitch_pp = pchip(w_ops,u_opsM(3,:));
        BP_fun = @(w) ppval(BPitch_pp,w);
        
        %% yo(t)*uo(t)
        % interpolate (Generator Speed*Generator Torque) into a continuous function
        GP = -u_opsM(2,:).*x_opsM(5,:);
        GP_pp = pchip(w_ops,GP);
        GP_fun = @(w) ppval(GP_pp,w);
        
    case 'Lagrange'
        % Setup the variables and matrices to calculate the lagrange term
        %%
        R1 = 1e-8;R2 = 1e+8;
        
        % -(yl+yo)(ul+uo) + Ru^2
        lx = 1;
        
        % Ru^2
        L(lx).left = 1;
        L(lx).right = 1;
        L(lx).matrix = diag([0,R1,R2]);
        lx = lx+1;
        
        % -yl*ul
        L(lx).left = 1;
        L(lx).right = 2;
        Lmat = zeros(3,101); Lmat(2,5) = 1;
        L(lx).matrix = -Lmat;
        lx = lx+1;
        
        % -ul*yo
        L(lx).left = 0;
        L(lx).right = 1;
        L2mat = cell(1,3);L2mat{2} = @(t) GSn_fun(W_fun(t));
        L(lx).matrix = L2mat;
        lx = lx+1;
        
        % -yl*uo
        L3mat = cell(1,101); L3mat{5} = @(t) GTn_fun(W_fun(t));
        L(lx).left = 0;
        L(lx).right = 2;
        L(lx).matrix = L3mat;
        lx = lx+1;
        
        
        TQmax = max(u_opsM(2,:));
        GSmax = 0.7913;
        SMat = zeros(101,1); SMat(5) = TQmax;
        
        % -yo*uo
        L(lx).left = 0;
        L(lx).right = 0;
        L(lx).matrix = {@(t) GP_fun(W_fun(t))};
        
        % linear constraint on power
        Z(1).linear(1).right = 1;
        Z(1).linear(1).matrix = [0,GSmax,0]';
        Z(1).linear(2).right = 2;
        Z(1).linear(2).matrix = SMat;
        Z(1).b = @(t)(15000 + 0.9655*TQmax*(GSmax + GSn_fun(W_fun(t))) + 0.9655*GSmax*(TQmax + GTn_fun(W_fun(t))));
        

    case 'IL_plot'
        % Innerloop plot code
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
        
end


