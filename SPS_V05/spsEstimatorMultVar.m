function [GammaHat,hyper] = spsEstimatorMultVar(X,Y,opt1,blk,process,hyper)

n=size(X,1);
d=size(X,2);
if size(Y,3)==1
    N=size(Y,2); r=1;
elseif size(Y,3)>1
    N=size(Y,3); r=size(Y,2);
end
% Estimating the precision matrix
[P,X,Y,~,Cn,Xn] = spsEstimator(X,Y,opt1,blk,process);
K=size(Cn,2);
C=Cn{1};
for k=2:K
    C=sparse(blkdiag(C,Cn{k}));
end
%% Estimating Gamma
display('Learning the between response covariance matrix ...');
GammaTilda=zeros(r,r,n);
for i=1:n-1
    beg=(i-1)*r+1;
    GammaTilda(:,:,i)=C(beg:beg+r-1,beg:beg+r-1);
end
GammaHat=mean(GammaTilda,3);
%% Estimating Spatial hyperparameters 
display('Learning hyper parameters ...');
hyper = covFuncEstimatorMultVar(Cn,X,hyper,GammaHat);
%% Closing
if strcmp(process.type,'parallel')
    matlabpool('close');
end

end