clear
clc
%% Initialization
n=100;
n0=100;
m=0;
N=10; %   No. of ensemble samples 
designRange=100;
phi=4; sigma2=1; tau2=0; mu=0;
% covFunc='Expiso'; p=[]; nu=[];
covFunc='Expiso'; p=[]; nu=[];
% covFunc='Materniso'; p=[]; nu=3/2;
%covFunc='PEiso'; p=4; nu=[];
%% Defining True Hyper Parameter
hyperTrue.covFunc=covFunc;
hyperTrue.p=p;
hyperTrue.nu=nu;
hyperTrue.param.val=[phi;sigma2;tau2];
if tau2==0
    hyperTrue.nugget='false';
else
    hyperTrue.nugget='true';
end
clear phi sigma2 tau2 covFunc p nu;
%% Building Design
X = designRange*lhsdesign(n,2);
X0 = designRange*lhsdesign(n0,2);
%% generate y
[~,~,~,~,Y,~]=GRFSim(X,hyperTrue,mu,N);
[~,~,~,~,Y0,~]=GRFSim(X0,hyperTrue,mu,1);
save('dataFileNPN.mat','X','Y','hyperTrue','mu','designRange','X0','Y0');