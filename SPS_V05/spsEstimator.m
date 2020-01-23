function [P,X,Y,identifier,Cn,Xn]= spsEstimator(X,Y,opt1,blk,process)
% written by Sam Davanloo-Tajbakhsh, 7/03/2015

n=size(Y,1);
N=size(Y,3);
r=size(Y,2);
muHatVec=mean(mean(Y,3),1)';
%% Blocking
identifier=[1:size(X,1)]';
if strcmp(blk.scheme,'SS') %Spatial Segmentation (SS)
    display('Blocking is based on Spatial Segmentation (SS)');
    for j=1:size(X,2)
        blkLims{j}=[blk.ss.min(j):blk.ss.range(j)/blk.ss.u(j):blk.ss.max(j)];
        blkLimLB{j}=blkLims{j}(1:size(blkLims{j},2)-1);
        blkLimUB{j}=blkLims{j}(2:size(blkLims{j},2));
    end
    lb=cartprod(blkLimLB{1:size(X,2)})';
    ub=cartprod(blkLimUB{1:size(X,2)})';
    k=1; h=1;
    while h<=blk.K
        LB=lb(:,h); UB=ub(:,h);
        condition='LB(1)<X(:,1)&X(:,1)<UB(1)';
        for j=2:size(X,2)
            condition=strcat(condition,'&','LB(',num2str(j),')<X(:,',num2str(j),')&X(:,',num2str(j),')<UB(',num2str(j),')');
        end
        [row,~]=find(eval(condition));
        if size(row,1)<=1 %<-------------------- could be set to a higher integer
            h=h+1;
            continue;
        end
        if r==1
            Xn{k}=X(row,:);
            Idt{k}=identifier(row,1);
            Yn{k}=Y(row,:);
            nn(k)=size(Xn{k},1);
            Dn{k}=squareform(pdist(Xn{k}));
            dMaxn(k)=max(max(Dn{k}));
            DTemp=Dn{k}+dMaxn(k)*eye(nn(k));
            minD=min(DTemp,[],2);
            Dpn{k}=Dn{k}+diag(minD);
%             lambda(k)=1/sqrt(nn(k));
            lambda(k)=500;
        elseif r>1
            nrow=size(row,1);
            Idt{k}=identifier(row,1);
            Xnew{k}=X(row,:);
            Xn{k}=[];
            for i=1:size(Xnew{k},1)
                Xn{k}=[Xn{k};repmat(Xnew{k}(i,:),r,1)];
            end
            Ynew{k}=Y(row,:,:);
            for l=1:N
                Yn{k}(:,l)=reshape(Y(row,:,l)',[nrow*r,1]);
            end
            nn(k)=size(Xn{k},1);
            Ds1=squareform(pdist(Xnew{k}));
            dMaxn(k)=max(max(Ds1));
            Ds2=Ds1+dMaxn(k)*eye(nrow);
            minD=min(Ds2,[],2);
            Dpn{k}=squareform(pdist(Xn{k}))+kron(diag(minD),ones(r,r));
            lambda(k)=1/sqrt(nn(k));
        end
        K=k;
        k=k+1;
        h=h+1;
    end
    blk.K=K;
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
        if r==1
            Idt{k}=identifier(row,1);
            Xn{k}=X(row,:);
            Yn{k}=Y(row,:);
            nn(k)=size(Xn{k},1);
            Dn{k}=squareform(pdist(Xn{k}));
            dMaxn(k)=max(max(Dn{k}));
            DTemp=Dn{k}+dMaxn(k)*eye(nn(k));
            minD=min(DTemp,[],2);
            Dpn{k}=Dn{k}+diag(minD);
            lambda(k)=1/sqrt(nn(k));
        elseif r>1
            Idt{k}=identifier(row,1);
            nrow=size(row,1);
            Xnew{k}=X(row,:);
            Xn{k}=[];
            for i=1:size(Xnew{k},1)
                Xn{k}=[Xn{k};repmat(Xnew{k}(i,:),r,1)];
            end
            Ynew{k}=Y(row,:,:);
            for l=1:N
                Yn{k}(:,l)=reshape(Y(row,:,l)',[nrow*r,1]);
            end
            nn(k)=size(Xn{k},1);
            Ds1=squareform(pdist(Xnew{k}));
            dMaxn(k)=max(max(Ds1));
            Ds2=Ds1+dMaxn(k)*eye(nrow);
            minD=min(Ds2,[],2);
            Dpn{k}=squareform(pdist(Xn{k}))+kron(diag(minD),ones(r,r));
            lambda(k)=1/sqrt(nn(k));
        end     
    end
end
%% Building X, D, Sigma and C0
display(sprintf('Max block size is %4d and Min block size is %4d',max(nn),min(nn)));
X=[]; Y=[]; identifier=[];
for k=1:K
    if r==1
        X=[X;Xn{k}];
        Y=[Y;Yn{k}];
    else
        X=[X;Xnew{k}];
        Y=[Y;Ynew{k}];
    end
    identifier=[identifier;Idt{k}];
end
%% First Optimization
if strcmp(process.type,'parallel')
    if (matlabpool('size')>0), matlabpool('close'), end;
    matlabpool('open',process.nCores);
    display(sprintf('Parallel processing using %d cores is selected',process.nCores));
    parfor k=1:K
        %muHat=mean(Yn{k}(:));
        muHat=kron(ones(n,N),muHatVec);
        Sn{k}=(Yn{k}-muHat)*(Yn{k}-muHat)'/N;
        Pn{k} = ADMM(Sn{k},Dpn{k}/dMaxn(k),lambda(k),opt1);
        Pn{k}=sparse(Pn{k});
        Cn{k}=Pn{k}\speye(size(Pn{k}));
    end
    matlabpool('close');
else % not parallel
    for k=1:K
        %muHat=mean(Yn{k}(:));
        muHat=kron(ones(n,N),muHatVec);
        Sn{k}=(Yn{k}-muHat)*(Yn{k}-muHat)'/N;
        Pn{k} = ADMM(Sn{k},Dpn{k}/dMaxn(k),lambda(k),opt1);
        Pn{k}=sparse(Pn{k});
        Cn{k}=Pn{k}\speye(size(Pn{k}));
    end
end
%% 
P=Pn{1};
for k=2:K
    P=blkdiag(P,Pn{k});
end

end