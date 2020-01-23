function [D,Sigma,Precision,X,Y,S]=GRFSim(X,hyperTrue,mu,N)
% Simulates and plots zero-mean GRF's with a powered exponential or matern covariance function
% Requires Statistics toolbox

n=size(X,1);

%% 2.
%Isotropic Case
D=squareform(pdist(X));
Sigma=covarianceFunc(hyperTrue,D);
Precision=inv(Sigma);

%% 3.
S=zeros(n,n);
Y=zeros(n,N);
for i=1:N
    y=mvnrnd(mu*ones(n,1),Sigma);
    y=y';
    Y(:,i)=y;
    S=S+y*y';  %Sample Covariance Matrix
end
S=S/N;
end