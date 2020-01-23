function hyper = covFuncEstimatorMultVar(Cn,X,hyper,Gamma)

K=size(Cn,2);
C=Cn{1};
for k=2:K
    C=sparse(blkdiag(C,Cn{k}));
end

if strcmp(hyper.covFunc,'AnisGeneral')
    lb=[hyper.param.LB.ell;hyper.param.LB.Delta(:)];
    ub=[hyper.param.UB.ell;hyper.param.UB.Delta(:)];
    theta0=ub.*rand(size(ub));
    options = optimset('Algorithm','interior-point');
    [thetaOpt,~]=fmincon(@(theta)stage2MultVarObjAnisGeneral(theta,C,Gamma,hyper,X),theta0,[],[],[],[],lb,ub,[],options);
    e=size(hyper.param.LB.ell,1);
    g=hyper.param.g;
    hyper.param.val.ell=thetaOpt(1:1+e-1);
    hyper.param.val.Delta=reshape(thetaOpt(1+e:1+e+(e*g)-1),[e,g]);
elseif strcmp(hyper.covFunc,'AnisDiag')
    lb=[hyper.param.LB.ell];
    ub=[hyper.param.UB.ell];
    theta0=ub.*rand(size(ub));
    options = optimset('Algorithm','interior-point');
    [thetaOpt,~]=fmincon(@(theta)stage2MultVarObjAnisDiag(theta,C,Gamma,hyper,X),theta0,[],[],[],[],lb,ub,[],options);
    e=size(hyper.param.LB.ell,1);
    hyper.param.val.ell=thetaOpt(1:1+e-1);   
else
    lb=hyper.param.LB(1);
    ub=hyper.param.UB(1);
    theta0=ub.*rand(size(ub));
    options = optimset('Algorithm','interior-point');
    [thetaOpt,~]=fmincon(@(theta)stage2MultVarObjIso(theta,C,Gamma,hyper,X),theta0,[],[],[],[],lb,ub,[],options);
    hyper.param.val(1)=thetaOpt;
    
end

end

function obj = stage2MultVarObjAnisGeneral(theta,C,Gamma,hyper,X)

e=size(hyper.param.LB.ell,1);
g=hyper.param.g;

hyper.param.val.ell=theta(1:1+e-1);
hyper.param.val.Delta=reshape(theta(1+e:1+e+(e*g)-1),[e,g]);

R=correlationFunc(hyper,X);
Sigma=kron(R,Gamma);

obj=(1/2)*norm(C-Sigma,'fro')^2;

end

function obj = stage2MultVarObjAnisDiag(theta,C,Gamma,hyper,X)

e=size(hyper.param.LB.ell,1);

hyper.param.val.ell=theta(1:1+e-1);

R=correlationFunc(hyper,X);
Sigma=kron(R,Gamma);

obj=(1/2)*norm(C-Sigma,'fro')^2;

end

function obj = stage2MultVarObjIso(theta,C,Gamma,hyper,X)

hyper.param.val(1)=theta;

R=correlationFunc(hyper,X);
Sigma=kron(R,Gamma);

obj=(1/2)*norm(C-Sigma,'fro')^2;

end

