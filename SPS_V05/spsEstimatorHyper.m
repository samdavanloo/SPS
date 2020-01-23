function [P,X,Y,hyper]= spsEstimatorHyper(X,Y,opt1,blk,process,hyper)
% written by Sam Davanloo-Tajbakhsh, 5/29/2014

if strcmp(process.type,'parallel')
    %if (matlabpool('size')>0), matlabpool('close'), end;
    matlabpool('open',process.nCores);
    display(sprintf('Parallel processing using %d cores is selected',process.nCores));
end

n=size(X,1);
if size(Y,3)==1
    N=size(Y,2); r=1;
elseif size(Y,3)>1
    N=size(Y,3); r=size(Y,2);
end
%% Blocking
if strcmp(blk.scheme,'SS') %Spatial Segmentation (SS)
    display('Blocking is based on Spatial Segmentation (SS)');
    for j=1:size(X,2)
        blkLims{j}=[blk.ss.min(j):blk.ss.range(j)/blk.ss.u(j):blk.ss.max(j)];
        blkLimLB{j}=blkLims{j}(1:size(blkLims{j},2)-1);
        blkLimUB{j}=blkLims{j}(2:size(blkLims{j},2));
    end
    lb=cartprod(blkLimLB{1:size(X,2)})';
    ub=cartprod(blkLimUB{1:size(X,2)})';
    
    k=0;
    for i=1:blk.K
        k=k+1;
        LB=lb(:,i); UB=ub(:,i);
        condition='LB(1)<X(:,1)&X(:,1)<UB(1)';
        for j=2:size(X,2)
            condition=strcat(condition,'&','LB(',num2str(j),')<X(:,',num2str(j),')&X(:,',num2str(j),')<UB(',num2str(j),')');
        end
        [row,~]=find(eval(condition));
        Xn{k}=X(row,:);
        Yn{k}=Y(row,:);
        nn(k)=size(Xn{k},1);
        Dn{k}=squareform(pdist(Xn{k}));
        dMaxn(k)=max(max(Dn{k}));
        DTemp=Dn{k}+dMaxn(k)*eye(nn(k));
        minD=min(DTemp,[],2);
        Dpn{k}=Dn{k}+diag(minD);
        lambda(k)=1/sqrt(nn(k));
    end

    K=blk.K;
else % Random Selection (RS)
    display('Blocking is based on Random Selection (RS)');
    K=blk.K;
    w=floor(n/K);
    nD=w*K;
    indc=1:nD;
    newIndc=indc(randperm(nD));
    newIndc=reshape(newIndc,w,K);
    for k=1:K
        if k < K
            row=newIndc(:,k);
        else
            row=[newIndc(:,k);[nD+1:n]'];
        end
        Xn{k}=X(row,:);
        Yn{k}=Y(row,:);
        nn(k)=size(Xn{k},1);
        Dn{k}=squareform(pdist(Xn{k}));
        dMaxn(k)=max(max(Dn{k}));
        DTemp=Dn{k}+dMaxn(k)*eye(nn(k));
        minD=min(DTemp,[],2);
        Dpn{k}=Dn{k}+diag(minD);
        lambda(k)=1/sqrt(nn(k));
    end
end
%% Building X, D, Sigma and C0
display(sprintf('Max block size is %4d and Min block size is %4d',max(nn),min(nn)));
X=[]; d=[]; Y=[];
for k=1:K
    X=[X;Xn{k}];
    Y=[Y;Yn{k}];
    DD=Dn{k};
    d=[d;DD(:)];
end
clear DD Xn;
%% First Optimization
%muHat=0; %<---------------------
if strcmp(process.type,'parallel')
    parfor k=1:K
        muHat=mean(Yn{k}(:));
        Sn{k}=(Yn{k}-muHat)*(Yn{k}-muHat)'/N;
        Pn{k} = ADMM(Sn{k},Dpn{k}/dMaxn(k),lambda(k),opt1);
        Pn{k}=sparse(Pn{k});
        Cn{k}=Pn{k}\speye(size(Pn{k}));
    end
else % not parallel
    for k=1:K
        muHat=mean(Yn{k}(:));
        Sn{k}=(Yn{k}-muHat)*(Yn{k}-muHat)'/N;
        Pn{k} = ADMM(Sn{k},Dpn{k}/dMaxn(k),lambda(k),opt1);
        Pn{k}=sparse(Pn{k});
        Cn{k}=Pn{k}\speye(size(Pn{k}));
    end
end
CC=Cn{1}; 
c=CC(:);
P=Pn{1};
for k=2:K
    CC=Cn{k};
    c=[c;CC(:)];
    P=blkdiag(P,Pn{k});
end
clear CC Cn Yn Pn Sn;
%% Second Optimization
display('Learning hyper parameters ...');
if strcmp(hyper.covFunc,'AnisGeneral')
    hyper = covFuncEstimatorAnisGeneral(c,hyper);
elseif strcmp(hyper.covFunc,'AnisDiag')
    hyper = covFuncEstimatorAnisDiag(c,hyper);
else
    thetaOpt = covFuncEstimator(c,d,hyper);
    if size(thetaOpt,1)==3
        hyper.param.val=thetaOpt;
    else
        hyper.param.val=[thetaOpt(1);thetaOpt(2);0];
    end
    clear c;
end
%% Closing
if strcmp(process.type,'parallel')
    matlabpool('close');
end

end