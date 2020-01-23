function f=stage2objAnisDiag(theta)

global c4s2 hyper4s2 X dm
n=size(X,1);
hyper4s2.param.val.var=theta(1);
hyper4s2.param.val.nugget=theta(2);
hyper4s2.param.val.ell=theta(3:3+dm-1);

Sigma = hyper4s2.param.val.var*correlationFunc(hyper4s2,X)+hyper4s2.param.val.nugget*eye(n);
f=norm(Sigma(:)-c4s2)^2;

end