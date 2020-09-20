function LinearModels = FASTLin_ProcessLinearModels(LinModelFile,FASTOutPath,ReduceModel,SaveFlag)
    % FAST parameters
    outPrefix = 'lin';
    outSuffix = '.outb';
    outFiles = dir(fullfile(FASTOutPath,[outPrefix,'*',outSuffix]));
    nLinCases = length(outFiles);
    if nLinCases <= 10
        numstring = '%01d';
    else
        numstring = '%02d';
    end
    
    % Initialize
    MBC = cell(1,nLinCases);
    matData = cell(1,nLinCases);
    WindSpeed = zeros(1,nLinCases);
    P = cell(1,nLinCases);
    
    % Temporary parameters
    iCase = 1;
    iAz = 1;
    
    % Read FST, Servo, and DISCON files
    FST = FAST2Matlab(fullfile(FASTOutPath,[outPrefix,'_',num2str(iCase-1,numstring),'.fst']));
    idx = find(contains(FST.Label,'ServoFile'));
    ServoFile = replace(FST.Val{idx},'"','');
    Servo = FAST2Matlab(fullfile(FASTOutPath,ServoFile));
    idx = find(contains(Servo.Label,'DLL_InFile'));
    DISCONFile = replace(Servo.Val{idx},'"','');
    DISCON = ROSCO2Matlab(fullfile(FASTOutPath,DISCONFile));
        
    % Obtain linearized model state, control, and output variables
    LinResult = ReadFASTLinear(fullfile(FASTOutPath,[outPrefix,'_',num2str(iCase-1,numstring),'.',num2str(iAz),'.lin']));
    nx = size(LinResult.x_op,1);
    nu = size(LinResult.u_op,1);
    ny = size(LinResult.y_op,1);
    xLabel = LinResult.x_desc;
    uLabel = LinResult.u_desc;
    yLabel = LinResult.y_desc;
    
    for iCase = 1:nLinCases
        % Read .lin files perform MBC3 transformation to obtain the linear
        % model independent to the azimuth angle of blades
        % Process .lin files
        LinFilesS{iCase} = dir(fullfile(FASTOutPath,[outPrefix,'_',num2str(iCase-1,numstring),'.*.lin']));
        if isempty(LinFilesS{1})
            error('WARNING: Didn''t find any linear files');
        end
        for iFile = 1:length(LinFilesS{iCase})
            LinFiles{iFile} = fullfile(FASTOutPath,LinFilesS{iCase}(iFile).name);
        end
        [MBC{iCase},matData{iCase}] = fx_mbc3(LinFiles);

        % Obtain FAST parameters and subfile parameters
        FSTName = fullfile(FASTOutPath,[outPrefix,'_',num2str(iCase-1,numstring),'.fst']);
        FP = FAST2Matlab(FSTName,2); % FAST parameters, 2 lines of header (v8)
        [IfWP, InflowFile] = GetFASTPar_Subfile(FP, 'InflowFile', FASTOutPath, FASTOutPath);
        [EdP, ElastoFile]  = GetFASTPar_Subfile(FP, 'EDFile', FASTOutPath, FASTOutPath);
        
        % State variable indices
        AzDesc = 'ED Variable speed generator DOF'; % Azimuth [rad]
        AzInd = find(contains(matData{iCase}.DescStates,AzDesc));
        indStates = 1:length(matData{iCase}.DescStates);
        indStates(AzInd) = []; % Remove azimuth state
        
        % Control variable indices
        WindDesc = 'IfW Extended input: horizontal wind speed';
        TqDesc = 'ED Generator torque';
        PitchDesc = 'ED Extended input: collective blade-pitch command';
        indInp.WindInd = find(contains(matData{iCase}.DescCntrlInpt,WindDesc));
        indInp.TqInd = find(contains(matData{iCase}.DescCntrlInpt,TqDesc));
        indInp.PitchInd = find(contains(matData{iCase}.DescCntrlInpt,PitchDesc));
        indInps = [indInp.WindInd, indInp.TqInd, indInp.PitchInd];
        indInpsLbl = {'IfW WindSpeedHor', 'ED GenTorq', 'ED BldPitchCommand'};
        
        % Output variable indices
        GenPwrDesc = 'SrvD GenPwr'; % Electrical generator power [kW]
        GenSpeedDesc = 'ED GenSpeed'; % Angular speed of high-speed shaft and generator [rpm]
        IPDeflDesc = 'ED IPDefl1'; % In-plane tip deflection of blade 1 (Same values for all blades) [m]
        NcIMURAxsDesc = 'ED NcIMURAxs'; % Nacelle IMU rotational acceleration in xs-axis [deg/s^2]
        NcIMURAysDesc = 'ED NcIMURAys'; % Nacelle IMU rotational acceleration in ys-axis [deg/s^2]
        NcIMURAzsDesc = 'ED NcIMURAzs'; % Nacelle IMU rotational acceleration in zs-axis [deg/s^2]
        NcIMUTAxsDesc = 'ED NcIMUTAxs'; % Nacelle IMU translational acceleration in xs-axis [deg/s^2]
        NcIMUTAysDesc = 'ED NcIMUTAys'; % Nacelle IMU translational acceleration in ys-axis [deg/s^2]
        NcIMUTAzsDesc = 'ED NcIMUTAzs'; % Nacelle IMU translational acceleration in zs-axis [deg/s^2]
        OoPDeflDesc = 'ED OoPDefl1'; % Out-of-plane tip deflection of blade 1 (Same values for all blades) [m]
        PlfmHeaveDesc = 'ED PtfmHeave'; % Platform vertical heave displacement [m]
        PlfmPitchDesc = 'ED PtfmPitch'; % Platform pitch tilt angular displacement [deg]
        PlfmRollDesc = 'ED PtfmRoll'; % Platform roll tilt angular displacement [deg]
        PlfmSurgeDesc = 'ED PtfmSurge'; % Platform horizontal surge displacement [m]
        PlfmSwayDesc = 'ED PtfmSway'; % Platform horizontal sway displacement [m]
        PlfmYawDesc = 'ED PtfmYaw'; % Platform yaw angular displacement [deg]
        RotThrustDesc = 'ED RotThrust'; % Rotor (low speed shaft) thrust force [kN]
        RotTqDesc = 'ED RotTorq'; % Rotor (low speed shaft) torque [kN-m]
        TTDspFADesc = 'ED TTDspFA'; % Tower top / yaw-bearing fore-aft deflection [m]
        TTDspSSDesc = 'ED TTDspSS'; % Tower top / yaw-bearing side-to-side deflection [m]
        TTDspTwstDesc = 'ED TTDspTwst'; % Tower top / yaw-bearing angular torsion deflection [deg]
        TwrBsFxtDesc = 'ED TwrBsFxt'; % Tower base fore-aft shear force [kN]
        TwrBsFytDesc = 'ED TwrBsFyt'; % Tower base side-to-side shear force [kN]
        TwrBsFztDesc = 'ED TwrBsFzt'; % Tower base axial force [kN]
        TwrBsMxtDesc = 'ED TwrBsMxt'; % Tower base roll (side-to-side) moment [kN-m]
        TwrBsMytDesc = 'ED TwrBsMyt'; % Tower base pitching (fore-aft) moment [kN-m]
        TwrBsMztDesc = 'ED TwrBsMzt'; % Tower base yaw moment [kN-m]
        YawBrFxpDesc = 'ED YawBrFxp'; % Tower top / yaw-bearing fore-aft shear force [kN]
        YawBrFypDesc = 'ED YawBrFyp'; % Tower top / yaw-bearing side-to-side shear force [kN]
        YawBrFzpDesc = 'ED YawBrFzp'; % Tower top / yaw-bearing axial force [kN]
        YawBrMxpDesc = 'ED YawBrMxp'; % Tower top / yaw-bearing roll moment [kN-m]
        YawBrMypDesc = 'ED YawBrMyp'; % Tower top / yaw-bearing pitch moment [kN-m]
        YawBrMzpDesc = 'ED YawBrMzp'; % Tower top / yaw-bearing yaw moment [kN-m]
        RtAeroCpDesc = 'AD RtAeroCp'; % Rotor aerodynamic power coefficient [-]
        RtAeroFxhDesc = 'AD RtAeroFxh'; % Rotor aerodynamic load in x-direction [N]
        RtAeroFyhDesc = 'AD RtAeroFyh'; % Rotor aerodynamic load in y-direction [N]
        RtAeroFzhDesc = 'AD RtAeroFzh'; % Rotor aerodynamic load in z-direction [N]
        RtAeroMxhDesc = 'AD RtAeroMxh'; % Rotor aerodynamic moment in x-direction [N-m]
        RtAeroMyhDesc = 'AD RtAeroMyh'; % Rotor aerodynamic moment in y-direction [N-m]
        RtAeroMzhDesc = 'AD RtAeroMzh'; % Rotor aerodynamic moment in z-direction [N-m]
        RtTSRDesc = 'AD RtTSR'; % Rotor tip-speed ratio [-]
        RtVAvgxhDesc = 'AD RtVAvgxh'; % Rotor-disk-averaged relative wind velocity in x-direction [m/s]
        Wave1ElevDesc = 'HD Wave1Elev'; % Wave motion [m]
        indOut.GenPwrInd = find(contains(matData{iCase}.DescOutput,GenPwrDesc));
        indOut.GenSpeedInd = find(contains(matData{iCase}.DescOutput,GenSpeedDesc));
        indOut.IPDeflInd = find(contains(matData{iCase}.DescOutput,IPDeflDesc));
        indOut.NcIMURAxsInd = find(contains(matData{iCase}.DescOutput,NcIMURAxsDesc));
        indOut.NcIMURAysInd = find(contains(matData{iCase}.DescOutput,NcIMURAysDesc));
        indOut.NcIMURAzsInd = find(contains(matData{iCase}.DescOutput,NcIMURAzsDesc));
        indOut.NcIMUTAxsInd = find(contains(matData{iCase}.DescOutput,NcIMUTAxsDesc));
        indOut.NcIMUTAysInd = find(contains(matData{iCase}.DescOutput,NcIMUTAysDesc));
        indOut.NcIMUTAzsInd = find(contains(matData{iCase}.DescOutput,NcIMUTAzsDesc));
        indOut.OoPDeflInd = find(contains(matData{iCase}.DescOutput,OoPDeflDesc));
        indOut.PlfmHeaveInd = find(contains(matData{iCase}.DescOutput,PlfmHeaveDesc));
        indOut.PlfmPitchInd = find(contains(matData{iCase}.DescOutput,PlfmPitchDesc));
        indOut.PlfmRollInd = find(contains(matData{iCase}.DescOutput,PlfmRollDesc));
        indOut.PlfmSurgeInd = find(contains(matData{iCase}.DescOutput,PlfmSurgeDesc));
        indOut.PlfmSwayInd = find(contains(matData{iCase}.DescOutput,PlfmSwayDesc));
        indOut.PlfmYawInd = find(contains(matData{iCase}.DescOutput,PlfmYawDesc));
        indOut.RotThrustInd = find(contains(matData{iCase}.DescOutput,RotThrustDesc));
        indOut.RotTqInd = find(contains(matData{iCase}.DescOutput,RotTqDesc));
        indOut.TTDspFAInd = find(contains(matData{iCase}.DescOutput,TTDspFADesc));
        indOut.TTDspSSInd = find(contains(matData{iCase}.DescOutput,TTDspSSDesc));
        indOut.TTDspTwstInd = find(contains(matData{iCase}.DescOutput,TTDspTwstDesc));
        indOut.TwrBsFxtInd = find(contains(matData{iCase}.DescOutput,TwrBsFxtDesc));
        indOut.TwrBsFytInd = find(contains(matData{iCase}.DescOutput,TwrBsFytDesc));
        indOut.TwrBsFztInd = find(contains(matData{iCase}.DescOutput,TwrBsFztDesc));
        indOut.TwrBsMxtInd = find(contains(matData{iCase}.DescOutput,TwrBsMxtDesc));
        indOut.TwrBsMytInd = find(contains(matData{iCase}.DescOutput,TwrBsMytDesc));
        indOut.TwrBsMztInd = find(contains(matData{iCase}.DescOutput,TwrBsMztDesc));
        indOut.YawBrFxpInd = find(contains(matData{iCase}.DescOutput,YawBrFxpDesc));
        indOut.YawBrFypInd = find(contains(matData{iCase}.DescOutput,YawBrFypDesc));
        indOut.YawBrFzpInd = find(contains(matData{iCase}.DescOutput,YawBrFzpDesc));
        indOut.YawBrMxpInd = find(contains(matData{iCase}.DescOutput,YawBrMxpDesc));
        indOut.YawBrMypInd = find(contains(matData{iCase}.DescOutput,YawBrMypDesc));
        indOut.YawBrMzpInd = find(contains(matData{iCase}.DescOutput,YawBrMzpDesc));
        indOut.RtAeroCpInd = find(contains(matData{iCase}.DescOutput,RtAeroCpDesc));
        indOut.RtAeroFxhInd = find(contains(matData{iCase}.DescOutput,RtAeroFxhDesc));
        indOut.RtAeroFyhInd = find(contains(matData{iCase}.DescOutput,RtAeroFyhDesc));
        indOut.RtAeroFzhInd = find(contains(matData{iCase}.DescOutput,RtAeroFzhDesc));
        indOut.RtAeroMxhInd = find(contains(matData{iCase}.DescOutput,RtAeroMxhDesc));
        indOut.RtAeroMyhInd = find(contains(matData{iCase}.DescOutput,RtAeroMyhDesc));
        indOut.RtAeroMzhInd = find(contains(matData{iCase}.DescOutput,RtAeroMzhDesc));
        indOut.RtTSRInd = find(contains(matData{iCase}.DescOutput,RtTSRDesc));
        indOut.RtVAvgxhInd = find(contains(matData{iCase}.DescOutput,RtVAvgxhDesc));
        indOut.Wave1ElevInd = find(contains(matData{iCase}.DescOutput,Wave1ElevDesc));
        indOuts = [indOut.GenPwrInd, indOut.GenSpeedInd, indOut.IPDeflInd, ...
            indOut.NcIMURAxsInd, indOut.NcIMURAysInd, indOut.NcIMURAzsInd, ...
            indOut.NcIMUTAxsInd, indOut.NcIMUTAysInd, indOut.NcIMUTAzsInd, indOut.OoPDeflInd, ...
            indOut.PlfmHeaveInd, indOut.PlfmPitchInd, indOut.PlfmRollInd, ...
            indOut.PlfmSurgeInd, indOut.PlfmSwayInd, indOut.PlfmYawInd, ...
            indOut.RotThrustInd, indOut.RotTqInd, ...
            indOut.TTDspFAInd, indOut.TTDspSSInd, indOut.TTDspTwstInd, ...
            indOut.TwrBsFxtInd, indOut.TwrBsFytInd, indOut.TwrBsFztInd, ...
            indOut.TwrBsMxtInd, indOut.TwrBsMytInd, indOut.TwrBsMztInd, ...
            indOut.YawBrFxpInd, indOut.YawBrFypInd, indOut.YawBrFzpInd, ...
            indOut.YawBrMxpInd, indOut.YawBrMypInd, indOut.YawBrMzpInd, indOut.RtAeroCpInd, ...
            indOut.RtAeroFxhInd, indOut.RtAeroFyhInd, indOut.RtAeroFzhInd, ...
            indOut.RtAeroMxhInd, indOut.RtAeroMyhInd, indOut.RtAeroMzhInd, ...
            indOut.RtTSRInd, indOut.RtVAvgxhInd, indOut.Wave1ElevInd];
        indOutsLbl = {GenPwrDesc, GenSpeedDesc, IPDeflDesc, ...
            NcIMURAxsDesc, NcIMURAysDesc, NcIMURAzsDesc, ...
            NcIMUTAxsDesc, NcIMUTAysDesc, NcIMUTAzsDesc, OoPDeflDesc, ...
            PlfmHeaveDesc, PlfmPitchDesc, PlfmRollDesc, ...
            PlfmSurgeDesc, PlfmSwayDesc, PlfmYawDesc, ...
            RotThrustDesc, RotTqDesc, ...
            TTDspFADesc, TTDspSSDesc, TTDspTwstDesc, ...
            TwrBsFxtDesc, TwrBsFytDesc, TwrBsFztDesc, ...
            TwrBsMxtDesc, TwrBsMytDesc, TwrBsMztDesc, ...
            YawBrFxpDesc, YawBrFypDesc, YawBrFzpDesc, ...
            YawBrMxpDesc, YawBrMypDesc, YawBrMzpDesc, RtAeroCpDesc, ...
            RtAeroFxhDesc, RtAeroFyhDesc, RtAeroFzhDesc, ...
            RtAeroMxhDesc, RtAeroMyhDesc, RtAeroMzhDesc, ...
            RtTSRDesc, RtVAvgxhDesc, Wave1ElevDesc};
        
        % Obtain operating point values
        SS_Op(iCase).xop = mean(matData{iCase}.xop,2);
        SS_Op(iCase).xdop = mean(matData{iCase}.xdop,2);
        SS_Op(iCase).uop = mean(matData{iCase}.uop,2);
        SS_Op(iCase).yop = mean(matData{iCase}.yop,2);
        
        % Reduce operating point vectors with indices of control and output variables
        SS_Ops(iCase).xop = SS_Op(iCase).xop(indStates',:);
        SS_Ops(iCase).xdop = SS_Op(iCase).xdop(indStates',:);
        SS_Ops(iCase).uop = SS_Op(iCase).uop(indInps',:);
        SS_Ops(iCase).yop = SS_Op(iCase).yop(indOuts',:);
        WindSpeed(iCase) = SS_Ops(iCase).uop(indInp.WindInd);
        
        % Form SS system
        P{iCase} = ss( ...
            MBC{iCase}.AvgA(indStates,indStates), ...
            MBC{iCase}.AvgB(indStates,indInps), ...
            MBC{iCase}.AvgC(indOuts,indStates), ...
            MBC{iCase}.AvgD(indOuts,indInps));
        % State Labels
        P{iCase}.StateName = MBC{iCase}.DescStates(indStates);
        % Input Labels
        P{iCase}.InputName = indInpsLbl;
        % Output Labels
        P{iCase}.OutputName = indOutsLbl;
        
        % Save original full size model
        P_full{iCase} = P{iCase};
        SS_Ops_full(iCase) = SS_Ops(iCase);
        
        % Model reduction
        if ReduceModel
            [P{iCase},UMat{iCase}] = minreal(P{iCase});
            % Operating point re-calculation for reduced model
            %funLSQ = @(x) (P{iCase}.C*x + P{iCase}.D*SS_Ops(iCase).uop - SS_Ops(iCase).yop)./(abs(SS_Ops(iCase).yop)+sqrt(eps));
            %xop = lsqnonlin(funLSQ,zeros(size(P{iCase}.A,1),1));
            xop = UMat{iCase}*SS_Ops_full(iCase).xop;
            xop = xop(1:size(P{iCase}.A,1),1);
            SS_Ops(iCase).xop = xop;
            SS_Ops(iCase).xdop = P{iCase}.A*xop + P{iCase}.B*SS_Ops(iCase).uop;
        end
    end
    
    % Return values
    LinearModels.P = P;
    LinearModels.P_full = P_full;
    LinearModels.UMat = UMat;
    LinearModels.SS_Ops = SS_Ops;
    LinearModels.SS_Ops_full = SS_Ops_full;
    LinearModels.WindSpeed = WindSpeed;
    LinearModels.FST = FST;
    LinearModels.Servo = Servo;
    LinearModels.DISCON = DISCON;
    
    % Save for future use
    if SaveFlag
        save(LinModelFile, 'P', 'P_full', 'UMat', 'SS_Ops', 'SS_Ops_full', 'WindSpeed', 'FST', 'Servo', 'DISCON');
    end
end
