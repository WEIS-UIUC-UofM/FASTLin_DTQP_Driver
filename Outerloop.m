% Outerloop.m
% Provides the constraints, wind inputs, linearized file to the innerloop
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
function Outerloop

clc; clear;

% pitch contraint
% Constraints = [3,4,5,6]; % list of pitch constraints
Constraints = 6;

% matfile corresponding to the PtfmMass
%CaseNames = {'pd_0.2_linear','pd_0.3_linear','pd_0.4_linear','pd_0.5_linear','pd_0.6_linear',...
%    'pd_0.7_linear','pd_0.8_linear','pd_0.9_linear',...
%    'pd_1.0_linear','pd_1.1_linear','pd_1.2_linear'};
CaseNames = 'pd_1.0_linear';

% Wind Files
%WindFiles = {'\iea15mw_0.mat','\iea15mw_1.mat','\iea15mw_2.mat','\iea15mw_3.mat','\iea15mw_4.mat','072720_183300.mat'};
WindFiles = '072720_183300.mat';

% calculate inner loop response
[T1,U,X,Pitch_rate,F] = Innerloop(CaseNames,Constraints,WindFiles);

end