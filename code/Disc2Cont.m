% Disc2Cont.m
% This script interpolates the discrete values of the state and control
% operating points into continuous functions dependednt on wind speed
% Specific functions of GeneratorSpeed, GeneratorTorque and BladePitch are
% required in the calculation of the objective term
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------

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
