function P = ADMM(S,D,mu,opt)

tolPrimal=opt.tol.primal;
tolDual=opt.tol.dual;
maxItr=opt.maxItr;
rho=opt.rho0;

n=size(S,1);
Z=eye(n,n);
U=zeros(n,n);
itr=0;
diagD=diag(D);
% m=2;
kappaI = 1.05;
% kappaD = 1.05;
s=Inf;
r=Inf;
normD=norm(D);

while (r>tolPrimal || s>tolDual) && itr<maxItr
    itr=itr+1;
           
    %Updating P
    [Q,L,~] = svd(rho*(Z-U)-S);
%     numRank=sum((diag(L)/rho)>1);
    lambda=diag(L);
    theta=(lambda+sqrt(lambda.^2+4*rho))/(2*rho);
    P=Q*diag(theta)*Q';
    
    %Updating Z
    Zold=Z;
    diagZ=max(diag(P+U)-mu*diagD/rho,0);
    Z = sign(P+U).*max(abs(P+U)-mu*D/rho,0);
    Z=Z.*(1-eye(size(Z)))+diag(diagZ);
    diffZ=Z-Zold;
    diffZmax(itr)=max(max(abs(diffZ)));
    
    %Updating rho
    PmZ=P-Zold;
    r=norm(PmZ(:))/normD;
    rr(itr)=r;
    s=rho*norm(diffZ)/normD;
    ss(itr)=s;
%     if r > m*s
%         rho=kappaI*rho;
%     elseif s > m*r
%         rho=rho/kappaD;
%     end
    rho = rho*kappaI;
    
    %Updating U
    Uold=U;
    U=(U+P-Z)/kappaI;
    diffU=U-Uold;
    diffUmax(itr)=max(max(abs(diffU)));
    
    % Obj
    logEigP=log(eig(P));
    logDet=sum(logEigP);
    trDZ=trace(D*abs(Z));
    objVal(itr)=trace(S'*P)-logDet+mu*trDZ;
    
    %Control
    if strcmp(opt.monitor,'on')
        if mod(itr,10)==0
            display(sprintf('itr:%3d, objVal: %1.1f, norm2PrResi: %1.2e, norm2DuResi: %1.2e',itr,objVal(itr),r,s));
        end
    end
    
end
display(sprintf('ADMM(Final)--> itr: %5d, objVal: %1.2f, norm2PrResi: %1.2e, norm2DuResi: %1.2e',itr,objVal(itr),r,s));

if strcmp(opt.monitor,'on')
    figure(1);
    subplot(3,1,1); plot([1:itr],log(rr(1:itr)),'b',[1:itr],log(ss(1:itr)),'r');
    h=legend('log(scaled primal residual)','log(scaled dual residual)');
    subplot(3,1,2); plot([1:itr],log(diffZmax(1:itr)),'b',[1:itr],log(diffUmax(1:itr)),'r');
    h=legend('log||Z_{k+1}-Z_k||_\infty','log||U_{k+1}-U_k||_\infty');
    subplot(3,1,3); plot([1:itr],objVal);
    h=legend('objective function');
    %figure(2); plot(diag(L)/rho,'r.'); ylabel('Singular Values');
end

P=Z;

end

