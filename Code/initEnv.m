%initEnv
% initialize environment
% Initialize environment
dbstop if error
% add the Framework folder
addpath(genpath(fullfile(pwd,'..','..','PolymerChainDynamics','Code')));
addpath(genpath(fullfile(pwd,'..','..','Utils')));