function thetaOpt = covFuncEstimator(C,D,hyper)

phiOptAux = @(x) phiOpt(x,C,D,hyper);

% Calculating the bounds
BS=hyper.param.UB(1);
fBS=phiOptAux(BS);
UB=BS+eps; LB=eps;

% Golden search
options=optimset('TolX',1e-4); %options=optimset('MaxFunEvals',20,'TolX',1e-2);
[theta1Opt,fOpt,exitflag,output] = fminbnd(phiOptAux,LB,UB,options);
[~,thetaOpt] = phiOpt(theta1Opt,C,D,hyper);
%display(thetaOpt);

end

function [fval,theta] = phiOpt(theta1,C,D,hyper)

if strcmp(hyper.nugget,'false')
    c=C;
    if ( strcmp(hyper.covFunc,'Expiso') || (strcmp(hyper.covFunc,'PEiso') && hyper.p==1) || (strcmp(hyper.covFunc,'Materniso') && hyper.nu==1/2) )
        aFunc = @(theta1,D) -exp(-D(:)./theta1);
    elseif ( strcmp(hyper.covFunc,'SEiso') || ( strcmp(hyper.covFunc,'PEiso') && hyper.p==2) )
        aFunc = @(theta1,D) -exp(-(D(:)./theta1).^2);
    elseif (strcmp(hyper.covFunc,'Materniso') && hyper.nu==3/2)
        aFunc = @(theta1,D) -(1+sqrt(3)*D(:)./theta1).*exp(-sqrt(3)*D(:)./theta1);
    elseif (strcmp(hyper.covFunc,'Materniso') && hyper.nu==5/2)
        aFunc = @(theta1,D) -(1+sqrt(5)*D(:)./theta1+5*D(:).^2/(3*theta1^2)).*exp(-sqrt(5)*D(:)./theta1);
    end
    
    a=aFunc(theta1,D);
    ac=a'*c;
    a2=norm(a)^2;
    
    if ac >= 0
        theta2=0;
    else
        theta2=-ac/a2;
    end
    
    fval=norm(c+a*theta2)^2;
    theta=[theta1;theta2];
else
    c=C;
    b=-(D(:)==0);
    if ( strcmp(hyper.covFunc,'Expiso') || (strcmp(hyper.covFunc,'PEiso') && hyper.p==1) || (strcmp(hyper.covFunc,'Materniso') && hyper.nu==1/2) )
        aFunc = @(theta1,D) -exp(-D(:)./theta1);
    elseif ( strcmp(hyper.covFunc,'SEiso') || ( strcmp(hyper.covFunc,'PEiso') && hyper.p==2) )
        aFunc = @(theta1,D) -exp(-(D(:)./theta1).^2);
    elseif (strcmp(hyper.covFunc,'Materniso') && hyper.nu==3/2)
        aFunc = @(theta1,D) -(1+sqrt(3)*D(:)./theta1).*exp(-sqrt(3)*D(:)./theta1);
    elseif (strcmp(hyper.covFunc,'Materniso') && hyper.nu==5/2)
        aFunc = @(theta1,D) -(1+sqrt(5)*D(:)./theta1+5*D(:).^2/(3*theta1^2)).*exp(-sqrt(5)*D(:)./theta1);
    end
    
    a=aFunc(theta1,D);
    bc=b'*c; ac=a'*c; ab=a'*b;
    b2=norm(b)^2; a2=norm(a)^2;
    theta2(1,1)=0;                            theta3(1,1)=-bc/b2;
    theta2(2,1)=-ac/a2;                       theta3(2,1)=0;
    theta2(3,1)=(ab*bc-ac*b2)/(a2*b2-(ab^2)); theta3(3,1)=(ab*ac-bc*a2)/(a2*b2-ab^2);
    theta2(4,1)=0;                            theta3(4,1)=0;
    
    if theta3(1,1) >= 0
        f(1,1) = norm(c+0+b*theta3(1,1))^2;
    else
        f(1,1) = Inf;
    end
    
    if theta2(2,1) >= 0
        f(2,1) = norm(c+a*theta2(2,1)+0)^2;
    else
        f(2,1) = Inf;
    end
    
    if ( theta2(3,1) >= 0 && theta3(3,1) >= 0 )
        f(3,1) = norm(c+a*theta2(3,1)+b*theta3(3,1))^2;
    else
        f(3,1) = Inf;
    end
    
    if ( theta2(4,1) == 0 && theta3(4,1) == 0 )
        f(4,1) = norm(c+0+0)^2;
    else
        f(3,1) = Inf;
    end
    
    [fval,indc]=min(f);
    theta=[theta1;theta2(indc);theta3(indc)];
end

end

