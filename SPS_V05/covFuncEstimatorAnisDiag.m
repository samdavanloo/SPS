function hyper = covFuncEstimatorAnisDiag(c,hyper)

lb=[hyper.param.LB.var;hyper.param.LB.nugget;hyper.param.LB.ell];
ub=[hyper.param.UB.var;hyper.param.UB.nugget;hyper.param.UB.ell];

theta0=ub.*rand(size(ub));

options = optimset('Algorithm','interior-point');

global c4s2 hyper4s2 dm
c4s2=c;
hyper4s2=hyper;

[thetaOpt,fVal]=fmincon(@stage2objAnisDiag,theta0,[],[],[],[],lb,ub,[],options);

hyper.param.val.var=thetaOpt(1);
hyper.param.val.nugget=thetaOpt(2);
hyper.param.val.ell=thetaOpt(3:3+dm-1);

end