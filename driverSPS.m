%%%% Written By: Sam Davanloo %%%%%
%%% This is the driver file for the SPS code V04. %%%
%%% For more information, please refer to http://arxiv.org/abs/1405.5576
close all
clear
clc
%% Add Path
addpath('C:\Users\sdt144\Dropbox\SamMATLAB\Sam\SPS_V04');  % add the path for the SPSV04 folder
%% Data
% Let:
% - d: dimension of the input space
% - n: number of data points
% - N: number of observed realizations of the data set
% Let X denote n*d matrix of input variables
% Let Y denote n*N matrix of response values

for j=1:d
    minX(j)=floor(min(X(:,j))); 
    maxX(j)=floor(max(X(:,j)))+1; 
    rangeX(j)=maxX(j)-minX(j); %finding the ranges of columns of X
end
%% Setting the segmentation scheme for big n
blk.scheme='SS'; % select either 'SS' (Spatial Segmentation) or 'RS' (Random Segmentation)   <---
if strcmp(blk.scheme,'SS') %% if SS-Blocking or No-Block is of interest
    uX=ones(d,1); % number of blocks along each coordinate; if no blocking is desired, set all equal to 1  <---
    blk.ss.min=minX;
    blk.ss.max=maxX;
    blk.ss.range=rangeX;
    blk.ss.u=uX; 
    blk.K=prod(uX); % total number of blocks for the SS blocking scheme
else % if RS-Blocking is of interest
    blk.K=1; % selected number of blocks for RS blocking scheme <---
end
%% Stage-I Optimization Parameters
opt1.method='ADMM'; % ADMM algorithm is used for stage-I optimization
opt1.monitor='off'; % 'off' or 'on'
opt1.tol.primal=1e-4; %primal feasibility threshold
opt1.tol.dual=1e-4; %dual feasibility threshold
opt1.maxItr=500; %number of ADMM iterations <---
opt1.rho0=(n/blk.K); %value of rho_0 in ADMM algorithm
%% Setting the covariance function to be used
hyper.covFunc='Materniso';%select the cov. func.:'SEiso' (squarede-exponential),'PEiso' (power-exponential),'Materniso'(matern),'Expiso' (exponential)<---
hyper.p=[]; % if power-exponential cov. function ('PEiso') is selected then hyper.p is the selected power, if not choose []; <---
hyper.nu=3/2; % if Matern cov. function ('Materniso') is selected then hyper.nu is the smoothing parameter, if not select []; <---
hyper.nugget='true'; %select from 'true' or 'false', if you want to add nugget to the covariance function <---
hyper.param.LB=[eps;eps;eps]; %lower bounds for [range;sigma2;nugget;nu;p]
dMax=norm(maxX-minX,2);
hyper.param.UB=[dMax;Inf;Inf]; %upper bound for [range;sigma2;nugget;nu;p]
%% Setting single vs. parallel processing
process.type='single'; % select 'parallel' if there is an access to multiple processing nodes or 'single' if not <---
if strcmp(process.type,'parallel')
    process.nCores=12; % define the number of available nodes for parallel processing <---
end
%% SPS estimation
[~,~,~,hyper] = spsEstimatorHyper(X,Y,opt1,blk,process,hyper);
        