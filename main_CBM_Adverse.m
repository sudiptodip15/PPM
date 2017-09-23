%% Projected Power Method (PPM) adapted for robust community detection. 
% The PPM algorithm is presented in the paper :
% ``The Projected Power Method: An Efficient Algorithm for Joint Alignment from Pairwise Differences'' by Y. Chen and E. J. Candes.
% **********************************************************************************************************************************
%  This modified code is for random corruption model with pairwise labels {+1,-1}. '+1' denotes two nodes are in same community 
%  and '-1' denotes they are in different community. '0' is for unobserved pairwise measurement.   
% **********************************************************************************************************************************

rng(0);

%% Set parameters
param.n             = 500;     % number of variables
param.m             = 10;      % number of states 
param.p_obs         = 0.7;     % observation ratio
param.maxIter       = 50;      % maximum number of projected power iterations
param.outFrac       = 0.6;     % Fraction of outliers
param.mu0           = 10;      % scale factor
param.pi0           = 0.75;    % non-corruption rate
param.round_flag    = 'no';    % indicate whether mu0 is infinity
param.max_poweriter = 20;      % maximum number of power iterations
param.eps_power     = 0.02;    % tolerance level of power methods
param.outlier_index = round(param.n*(1-param.outFrac)); % Compute outlier index
param.adv_model = 'rnd';
alg_num = 1;
param.max_robust_iter = 10;

% Noise density for the random corruption model
P0          = (1-param.pi0)* ones( param.m, 1) / param.m; 
P0(1)       = P0(1) + param.pi0;
param.P0    = P0; 

%% Generate a synthetic example

Problem=ProbGenAdverse(param);                     
fprintf('Running Robust Alg -->\n');     % Robust Cluster Assignment motivated by robust regression
RobustPPM_Alg(Problem,param); 
   

