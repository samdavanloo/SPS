function f=stage2objAnisGeneral(theta)

global c4s2 hyper4s2 X dm dk
n=size(X,1);
hyper4s2.param.val.var=theta(1);
hyper4s2.param.val.nugget=theta(2);
hyper4s2.param.val.ell=theta(3:3+dm-1);
hyper4s2.param.val.Delta=reshape(theta(3+dm:3+dm+(dm*dk)-1),[dm,dk]);

Sigma = hyper4s2.param.val.var*correlationFunc(hyper4s2,X)+hyper4s2.param.val.nugget*eye(n);
f=norm(Sigma(:)-c4s2)^2;

end