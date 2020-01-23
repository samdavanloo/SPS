function Sigma = covarianceFunc(hyper,X)

n=size(X,1);
covFunc=hyper.covFunc;
if strcmp(covFunc,'AnisDiag') || strcmp(covFunc,'AnisGeneral')
    sigma2=hyper.param.val.var;
else
    p=hyper.p;
    nu=hyper.nu;
    phi=hyper.param.val(1);
    sigma2=hyper.param.val(2);
end


if (strcmp(covFunc,'Expiso') || (strcmp(covFunc,'PEiso') && p==1) || (strcmp(covFunc,'Materniso') && nu==1/2))
    D=squareform(pdist(X));
    R=exp(-D/phi);
elseif (strcmp(covFunc,'SEiso') || (strcmp(covFunc,'PEiso') && p==2))
    D=squareform(pdist(X));
    R=exp(-(D/phi).^2);
elseif strcmp(covFunc,'PEiso')
    D=squareform(pdist(X));
    R=exp(-(D/phi).^p);
elseif (strcmp(covFunc,'Materniso') && (nu==3/2))
    D=squareform(pdist(X));
    R=(1+sqrt(3)*D/phi).*exp(-sqrt(3)*D/phi);
elseif (strcmp(covFunc,'Materniso') && (nu==5/2))
    D=squareform(pdist(X));
    R=(1+sqrt(5)*D/phi+5*D.^2/(3*phi^2)).*exp(-sqrt(5)*D/phi);
elseif strcmp(covFunc,'Materniso')
    D=squareform(pdist(X));
    Temp1=D/phi;
    Temp2=(1/(2^(nu-1)*gamma(nu)))*(Temp1.^nu);
    R=Temp2.*besselk(nu,Temp1);
    R(isnan(R))=1;
elseif strcmp(covFunc,'PEanis')
    D=squareform(pdist(X));
    sigma2=hyper.param.val(1);
    L=size(hyper.param.val,1)-2;
    phi=hyper.param.val(3:L+2);
    Temp=zeros(n,n);
    for j=1:L
        Temp=Temp+(D(:,:,j)/phi(j)).^p;
    end
    R=exp(-Temp);
elseif strcmp(covFunc,'AnisDiag')
    M=diag(hyper.param.val.ell.^(-1));
    D=mahal(X,M);
    R =exp(-D);
elseif strcmp(covFunc,'AnisGeneral')
    M=diag(hyper.param.val.ell.^(-1))+hyper.param.val.Delta*hyper.param.val.Delta';
    %D=squareform(pdist(X,'mahalanobis',inv(hyper.param.val.M)));
    %R =exp((D.^2)/(-2));
    D=mahal(X,M);
    R =exp(-D);
else
    error('The covariance function is undefined');
end

Sigma=sigma2*R;

end