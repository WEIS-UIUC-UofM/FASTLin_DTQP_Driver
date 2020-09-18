function FASTLin_DTQP_Test01

    clc;
    mname = mfilename('fullpath');
    [mpath, mname] = fileparts(mname);
    savepath = fullfile(mpath,'SaveData',mname);
    
    % Read/Create Linear Models
    LinModelFile = fullfile(savepath,'LinearModels.mat');
    if (exist(LinModelFile,'file') == 2)
        LinearModels = load(LinModelFile);
    else
        FASTOutPath = '/home/weis/ylee196/SaveData/TrimTest/LinearTwrPit_Tol1en5';
        ReduceModel = 1;
        SaveFlag = 1;
        LinearModels = FASTLin_ProcessLinearModels(LinModelFile,FASTOutPath,ReduceModel,SaveFlag);
    end
    
    % Read Wind Disturbance Profile
    WindFile = fullfile(savepath,'072720_183300.mat');
    if (exist(WindFile,'file') == 2)
        % Reading wind profile
        Wind_o = load(WindFile);
        Wind_o.NL_startTime = 450;
        Wind_o.Lin_TMax = 150;
        lin_inds = Wind_o.Chan.tt >= Wind_o.NL_startTime ...
            & Wind_o.Chan.tt <= Wind_o.NL_startTime + Wind_o.Lin_TMax;
        % Wind profile preparation
        Wind_data = Wind_o.Chan.RtVAvgxh(lin_inds);
        Wind_mean = mean(Wind_data);
        Wind_time = Wind_o.Chan.tt(lin_inds);
        Wind_time = Wind_time - Wind_time(1);
        Wind_sampling_dt = 0.5;
        Wind_sampling_nt = ceil(Wind_time(end)/Wind_sampling_dt) + 1;
        Wind_data_smooth = smoothdata(Wind_data,'sgolay',Wind_sampling_nt);
        Wind_profile = @(t) interp1(Wind_time,Wind_data_smooth,t);
    else
        error('Wind disturbance profile not found');
    end
    
    % Creating interpolated linear model
    AA = zeros(length(LinearModels.P{1}.A),length(LinearModels.P{1}.A),length(LinearModels.P));
    BB = zeros(size(LinearModels.P{1}.B,1),size(LinearModels.P{1}.B,2),length(LinearModels.P));
    CC = zeros(size(LinearModels.P{1}.C,1),size(LinearModels.P{1}.C,2),length(LinearModels.P));
    DD = zeros(size(LinearModels.P{1}.D,1),size(LinearModels.P{1}.D,2),length(LinearModels.P));
    u_opsM  = zeros(size(LinearModels.SS_Ops(1).uop,1),length(LinearModels.P));
    y_opsM  = zeros(size(LinearModels.SS_Ops(1).yop,1),length(LinearModels.P));
    for iSys = 1:length(LinearModels.P)
        s = size(LinearModels.P{iSys}.A);
        nStates(iSys) = s(1);
        uh_ops(iSys) = LinearModels.WindSpeed(iSys);
        u_opsM(:,iSys) = LinearModels.SS_Ops(iSys).uop;
        y_opsM(:,iSys) = LinearModels.SS_Ops(iSys).yop;

        AA(:,:,iSys)    = LinearModels.P{iSys}.A;
        BB(:,:,iSys)    = LinearModels.P{iSys}.B;
        CC(:,:,iSys)    = LinearModels.P{iSys}.C;
        DD(:,:,iSys)    = LinearModels.P{iSys}.D;
    end
    AA_p = permute(AA,[3,1,2]);
    BB_p = permute(BB,[3,1,2]);
    CC_p = permute(CC,[3,1,2]);
    DD_p = permute(DD,[3,1,2]);
    
    % Interpolate System Matrices & Operating Points
    A_op = squeeze(interp1(uh_ops,AA_p,Wind_mean));
    B_op = squeeze(interp1(uh_ops,BB_p,Wind_mean));
    C_op = squeeze(interp1(uh_ops,CC_p,Wind_mean));
    D_op = squeeze(interp1(uh_ops,DD_p,Wind_mean));
    P_op = ss(A_op,B_op,C_op,D_op);
    P_op.InputName = LinearModels.P{1}.InputName;
    P_op.OutputName = LinearModels.P{1}.OutputName;
    u_op = interp1(uh_ops',u_opsM',Wind_mean);
    y_op = interp1(uh_ops',y_opsM',Wind_mean);
    %x_op cannot be computed, since the model is reduced
    
    % DT-QP
    opts.general.mname = mname;
    opts.general.mpath = mpath;
    opts.general.name = mname;
    opts.general.path = savepath;
    opts.general.plotflag = 1;  % 0:do not plot, 1:plot
    opts.general.saveflag = 1;  % 0:do not save, 1:save
    opts.general.displevel = 2; % 0:none, 1:minimal, 2:verbose
    
    opts.dt.defects = 'TR';     % ZO, EF, Huen, ModEF, TR, HS, RK4, PS
    opts.dt.quadrature = 'CTR'; % CEF, CTR, CQHS, G, CC
    opts.dt.mesh = 'ED';        % ED, LGL, CGL, USER
    opts.dt.nt = 1000;
    opts.dt.meshr.method = 'none'; % none, richardson-doubling, ss-betts
    
    opts.method.reordervariables = 0;
    
    opts.solver.function = 'quadprog';
    opts.solver.tolerance = 1e-15;
    opts.solver.maxiter = 1000;
    opts.solver.display = 'none';
    
    % Objective function = maximize energy (integral of generator power)
    % ( Later if needed we may move to (1/2)*rho*A_sweep*(V_relative)^3,
    %   where V_relative = V_wind - V_towertop )
    idx = find(contains(P_op.OutputName,'SrvD GenPwr'));
    L(1).left = 0;
    L(1).right = 2;
    L(1).matrix = -P_op.C(idx,:);
    L(2).left = 0;
    L(2).right = 1;
    L(2).matrix = -P_op.D(idx,:);
    
    % Constraint function 1: rated power
    Z(1).linear(1).right = 2;
    Z(1).linear(1).matrix = -L(1).matrix;
    Z(1).linear(2).right = 1;
    Z(1).linear(2).matrix = -L(2).matrix;
    Z(1).b = 15000; % 15MW rated power
    
    % Upper bounds for states
    UB(1).right = 2;
    UB(1).matrix = 100*ones(size(A_op,1),1);
    
    % Upper bounds for controls
    idx1 = find(contains(LinearModels.DISCON.Label,'VS_RtTq'));
    idx2 = find(contains(LinearModels.DISCON.Label,'SD_MaxPit'));
    UB(2).right = 1;
    UB(2).matrix = {@(t)Wind_profile(t);
                    LinearModels.DISCON.Val{idx1};
                    LinearModels.DISCON.Val{idx2}};
    
    % Upper bounds for states at initial time
    % YHL: If we do not have this initial condition, several problems occur
    % including constraint violation, large discontinuity, etc.
    %UB(3).right = 4;
    %UB(3).matrix = zeros(size(A_op,1),1);
    
    % Lower bounds for states
    LB(1).right = 2;
    LB(1).matrix = -100*ones(size(A_op,1),1);
    
    % Lower bounds for controls
    LB(2).right = 1;
    LB(2).matrix = {@(t)Wind_profile(t);
                    0;
                    0};
    
    % Lower bounds for states at initial time
    % YHL: If we do not have this initial condition, several problems occur
    % including constraint violation, large discontinuity, etc.
    %LB(3).right = 4;
    %LB(3).matrix = zeros(size(A_op,1),1);
    
    setup.A = P_op.A;
    setup.B = P_op.B;
    setup.C = P_op.C;
    setup.D = P_op.D;
    setup.L = L;
    setup.Z = Z;
    setup.UB = UB;
    setup.LB = LB;
    setup.tf = Wind_time(end);
%    setup.p = p;
    
    [T1,U1,Y1,P1,F1,in,opts] = DTQP_solve(setup,opts);
    tic;
    [R2,T2,Y2] = lsim(P_op,U1,T1);
    fprintf('---------------------------------------------------\n');
    fprintf('simulation time for the optimal solution: %8.6f\n',toc);
    
    % Plot
    close all;
    fg1 = figure('Color',[1 1 1]);
    fg1.Position = [50 100 1500 800];
    % Plot Control
    subplot(6,3,1);
        idx = find(contains(P_op.InputName,'IfW WindSpeedHor'));
        plot(T1,U1(:,idx));
        title('Horizontal Wind Speed [m/s]');
    subplot(6,3,4);
        idx = find(contains(P_op.InputName,'ED GenTorq'));
        plot(T1,U1(:,idx)/1000);
        title('Generator Torque [kN-m]');
    subplot(6,3,7);
        idx = find(contains(P_op.InputName,'ED BldPitchCommand'));
        plot(T1,U1(:,idx)/pi*180);
        title('Blade Pitching Control [deg]');
        
    subplot(6,3,10);
        idx1 = find(contains(P_op.OutputName,'ED GenSpeed'));
        idx2 = find(contains(P_op.OutputName,'ED RotTorq'));
        plot(T2,R2(:,idx1).*R2(:,idx2));
        title('Gen speed * Rot Torque');
    % Plot Output
    subplot(6,3,2);
        idx = find(contains(P_op.OutputName,'SrvD GenPwr'));
        plot(T2,R2(:,idx));
        title('Generator Power [kW]');
        ylim([0 18000]);
        yticks([0 3000 6000 9000 12000 15000 18000]);
    subplot(6,3,5);
        idx = find(contains(P_op.OutputName,'ED GenSpeed'));
        plot(T2,R2(:,idx));
        title('Generator Speed [rpm]');
        ylim([0 8]);
        yticks([0 2 4 6 8]);
    subplot(6,3,8);
        idx = find(contains(P_op.OutputName,'ED IPDefl1'));
        tmp1 = plot(T2,R2(:,idx)); hold on;
        idx = find(contains(P_op.OutputName,'ED OoPDefl1'));
        tmp2 = plot(T2,R2(:,idx));
        legend([tmp1,tmp2],{'in-plane','out-of-plane'});
        title('Blade Deflection [m]');
    subplot(6,3,11);
        idx = find(contains(P_op.OutputName,'ED NcIMURAxs'));
        tmp1 = plot(T2,R2(:,idx)); hold on;
        idx = find(contains(P_op.OutputName,'ED NcIMURAys'));
        tmp2 = plot(T2,R2(:,idx));
        idx = find(contains(P_op.OutputName,'ED NcIMURAzs'));
        tmp3 = plot(T2,R2(:,idx));
        legend([tmp1,tmp2,tmp3],{'x-axis','y-axis','z-axis'});
        title('Nacelle IMU Angular Acceleration [deg/s^2]');
    subplot(6,3,14);
        idx = find(contains(P_op.OutputName,'ED NcIMUTAxs'));
        tmp1 = plot(T2,R2(:,idx)); hold on;
        idx = find(contains(P_op.OutputName,'ED NcIMUTAys'));
        tmp2 = plot(T2,R2(:,idx));
        idx = find(contains(P_op.OutputName,'ED NcIMUTAzs'));
        tmp3 = plot(T2,R2(:,idx));
        legend([tmp1,tmp2,tmp3],{'x-axis','y-axis','z-axis'});
        title('Nacelle IMU Translational Acceleration [m/s^2]');
    subplot(6,3,17);
        idx = find(contains(P_op.OutputName,'ED PtfmSurge'));
        tmp1 = plot(T2,R2(:,idx)); hold on;
        idx = find(contains(P_op.OutputName,'ED PtfmSway'));
        tmp2 = plot(T2,R2(:,idx));
        idx = find(contains(P_op.OutputName,'ED PtfmHeave'));
        tmp3 = plot(T2,R2(:,idx));
        legend([tmp1,tmp2,tmp3],{'surge','sway','heave'});
        title('Platform Translation [m]');
    subplot(6,3,3);
        idx = find(contains(P_op.OutputName,'ED PtfmRoll'));
        tmp1 = plot(T2,R2(:,idx)); hold on;
        idx = find(contains(P_op.OutputName,'ED PtfmPitch'));
        tmp2 = plot(T2,R2(:,idx));
        idx = find(contains(P_op.OutputName,'ED PtfmYaw'));
        tmp3 = plot(T2,R2(:,idx));
        legend([tmp1,tmp2,tmp3],{'roll','pitch','yaw'});
        title('Platform Tilt [deg]');
    subplot(6,3,6);
        idx = find(contains(P_op.OutputName,'ED RotThrust'));
        plot(T2,R2(:,idx));
        title('Rotor Thrust [kN]');
    subplot(6,3,9);
        idx = find(contains(P_op.OutputName,'ED RotTorq'));
        plot(T2,R2(:,idx));
        title('Rotor Torque [kN-m]');
    subplot(6,3,12);
        idx = find(contains(P_op.OutputName,'AD RtAeroCp'));
        plot(T2,R2(:,idx));
        title('Aerodynamic Power Coefficient [-]');
    subplot(6,3,15);
        idx = find(contains(P_op.OutputName,'AD RtTSR'));
        plot(T2,R2(:,idx));
        title('Rotor Tip Speed Ratio [-]');
    subplot(6,3,18);
        idx = find(contains(P_op.OutputName,'HD Wave1Elev'));
        plot(T2,R2(:,idx));
        title('Hydrodynamic Wave Elevation [m]');
end