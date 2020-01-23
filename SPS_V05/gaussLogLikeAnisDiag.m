function l = gaussLogLikeAnisDiag(theta)

global Y hyperMLE X dm

n=size(Y,1);

hyperMLE.param.val.var=theta(1);
hyperMLE.param.val.nugget=theta(2);
hyperMLE.param.val.ell=theta(3:3+dm-1);

Sigma=covarianceFunc(hyperMLE,X)+hyperMLE.param.val.nugget*eye(n);
if sum(eig(Sigma)>0)<n
    notPD=1;
end
e=0;
P=inv(Sigma+e*eye(n));
N=size(Y,2);
t=0;
for j=1:N
    t=t+Y(:,j)'*P*Y(:,j);
end
C=chol(Sigma+e*eye(n));

%d=prod(diag(C).^2);
%l=log(d)+t/N;

d=log(diag(C).^2);
l=sum(d)+t/N;

end