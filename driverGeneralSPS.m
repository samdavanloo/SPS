%%%% Written By: Sam Davanloo %%%%%
%%% This is the driver file for the General SPS code Version 05 (SPSV05)%%%
%%% For more information, please refer to ---------------------
close all
clear
clc
%% Add Path
% Add the folder 'SPSV05' to the path
%addpath('C:\Users\sdt144\Dropbox\SamMATLAB\Sam\SPS_V05');
addpath('/Users/Sam/Dropbox/SamMATLAB/Sam/SPS_V05');
%% Data
% X denotes n*d matrix of input variables
% Y denotes n*N*1 univariate or n*r*N r-variate response

% where:
% - d: dimension of the input space
% - n: number of data points
% - N: number of observed realizations of the data set
% - r: number of response variables

% Include a line that reads the data {X,Y} <----------------
X=randn(30,2);
Y=randn(30,2,10);

n=size(X,1);
d=size(X,2);
if size(Y,3)==1
    N=size(Y,2); r=1;
elseif size(Y,3)>1
    N=size(Y,3); r=size(Y,2);
end

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
opt1.monitor='off'; % 'off' or 'on' <---
opt1.tol.primal=1e-4; %primal feasibility threshold
opt1.tol.dual=1e-4; %dual feasibility threshold
opt1.maxItr=500; %number of ADMM iterations <---
opt1.rho0=(n/blk.K); %value of rho_0 in ADMM algorithm
%% Setting the covariance function to be used
hyper.covFunc='AnisGeneral'; % select from 'SEiso'(squarede-exponential),... <---
% 'PEiso'(power-exponential),'Materniso'(Matern),'Expiso'(exponential),...
% 'AnisDiag'(Diagonal Anisotropic), or 'AnisGeneral'(General Anisotropic)
if ( strcmp(hyper.covFunc,'AnisDiag') || strcmp(hyper.covFunc,'AnisGeneral') )
    % for more information see the "Gaussian Process for Machine Leaning by Rasmussen page 106.
    % Upper Bounds could be different for different parameters. 
    hyper.param.LB.var=eps;
    hyper.param.UB.var=100;
    if strcmp(hyper.covFunc,'AnisDiag')
        hyper.param.LB.ell=eps*ones(d,1);
        hyper.param.UB.ell=100*ones(d,1);
    elseif strcmp(hyper.covFunc,'AnisGeneral')
        hyper.param.LB.ell=eps*ones(d,1);
        g=1; % g is a positive integer number less than d (g<d)
        hyper.param.g=g;
        hyper.param.LB.ell=eps*ones(d,1);
        hyper.param.UB.ell=100*ones(d,1);
        hyper.param.LB.Delta=eps*ones(d,g);
        hyper.param.UB.Delta=1*ones(d,g); % should be adjusted carefully
    end
    hyper.param.LB.nugget=eps; % only used for univariate case
    hyper.param.UB.nugget=100; % only used for univariate case
else
    hyper.p=[]; % if power-exponential cov. function ('PEiso') is selected then hyper.p is the selected power, if not choose []; <---
    hyper.nu=3/2; % if Matern cov. function ('Materniso') is selected then hyper.nu is the smoothing parameter, if not select []; <---
    hyper.nugget='true'; %select from 'true' or 'false', if you want to add nugget to the covariance function <---
    hyper.param.LB=[eps;eps;eps]; %lower bounds for [range;sigma2;nugget]
    dMax=norm(maxX-minX,2);
    hyper.param.UB=[dMax;Inf;Inf]; %upper bound for [range;sigma2;nugget]
end
%% Setting single vs. parallel processing
process.type='single'; % select 'parallel' if there is an access to multiple processing nodes or 'single' if not <---
if strcmp(process.type,'parallel')
    process.nCores=12; % define the number of available nodes for parallel processing <---
end
%% SPS estimation
if r==1 % univariate response
    [~,~,~,hyper] = spsEstimatorHyper(X,Y,opt1,blk,process,hyper); % Estimation
    if ~( strcmp(hyper.covFunc,'AnisDiag') || strcmp(hyper.covFunc,'AnisGeneral') )
        rangeHat=hyper.param.val(1);
        varianceHat=hyper.param.val(2);
        nuggetHat=hyper.param.val(2);
        display(sprintf('range: %0.2f, variance: %0.2f, nugget: %0.2f',rangeHat,varianceHat,nuggetHat));
    elseif strcmp(hyper.covFunc,'AnisDiag')
        Mhat=diag(hyper.param.val.ell.^(-1)); display(Mhat);
        varianceHat=hyper.param.val.var; display(varianceHat);
        nuggetHat=hyper.param.val.nugget; display(nuggetHat);
    elseif strcmp(hyper.covFunc,'AnisGeneral')
        Mhat=diag(hyper.param.val.ell.^(-1))+hyper.param.val.Delta*hyper.param.val.Delta'; display(Mhat);
        varianceHat=hyper.param.val.var; display(varianceHat);
        nuggetHat=hyper.param.val.nugget; display(nuggetHat);
    end
else
    [GammaHat,hyper] = spsEstimatorMultVar(X,Y,opt1,blk,process,hyper); % Estimation
    display(GammaHat);
    if ~( strcmp(hyper.covFunc,'AnisDiag') || strcmp(hyper.covFunc,'AnisGeneral') )
        rangeHat=hyper.param.val(1); display(rangeHat);
    elseif strcmp(hyper.covFunc,'AnisDiag')
        Mhat=diag(hyper.param.val.ell.^(-1)); display(Mhat);
    elseif strcmp(hyper.covFunc,'AnisGeneral')
        Mhat=diag(hyper.param.val.ell.^(-1))+hyper.param.val.Delta*hyper.param.val.Delta'; display(Mhat);
    end
    
end
