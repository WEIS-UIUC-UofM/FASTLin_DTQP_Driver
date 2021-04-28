% Lagrange.m
% Setup the variables and matrices to calculate the lagrange term
%--------------------------------------------------------------------------
% 
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
% 
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


