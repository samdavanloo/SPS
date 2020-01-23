function l = gaussLogLike(theta)

global D Y hyperMLE;

hyperLocal=hyperMLE;
hyperLocal.param.val=theta;
n=size(D,1);
Sigma=covarianceFunc(hyperLocal,D)+hyperLocal.param.val(3)*eye(n);
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