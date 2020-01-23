function g = powerTrans(y,mu,sigma)

g0 = @(t) sign(t).*(abs(t)).^3;
denom3 = @(t,mu,sigma) ((g0(t-mu)).^2).*normpdf((t-mu)/sigma);
denom1 = @(t) denom3(t,mu,sigma);

denomVal = integral(denom1,-Inf,Inf);
g = sigma*(g0(y-mu)/sqrt(denomVal))+mu;

end