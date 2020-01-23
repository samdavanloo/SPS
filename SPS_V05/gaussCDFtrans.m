function g = gaussCDFtrans(y,mu,sigma)

global int11Val

g0 = @(t) normcdf((t-0.05)/0.4);
integrand13 = @(t,mu,sigma) ((g0(t)).^2).*normpdf((t-mu)/sigma);
integrand11 = @(t) integrand13(t,mu,sigma);
integrand24 = @(t,mu,sigma,int11Val) ((g0(t)-int11Val).^2).*normpdf((t-mu)/sigma);
integrand21 = @(t) integrand24(t,mu,sigma,int11Val);


int11Val = integral(integrand11,-Inf,Inf);
int21Val = integral(integrand21,-Inf,Inf);
g = sigma*((g0(y)-int11Val)/sqrt(int21Val))+mu;

end


% function g = gaussCDFtrans(z,M,S)
% 
% global mu sigma int1
% 
% mu=M;
% sigma=S;
% int1 = integral(@intgrndNum,-Inf,Inf);
% int2 = integral(@intgrndDenom,-Inf,Inf);
% g = sigma*((g0(z)-int1)/sqrt(int2))+mu;
% 
% end
% 
% function g_0 = g0(t)
% g_0 = normcdf((t-0.05)/0.4);
% end
% 
% function intgrndNumVal = intgrndNum(t)
% global mu sigma
% intgrndNumVal = ((g0(t)).^2).*normpdf((t-mu)/sigma);
% end
% 
% function intgrndDenomVal = intgrndDenom(y)
% global mu sigma int1
% intgrndDenomVal = (g0(y)-int1).^2.*normpdf((y-mu)/sigma);
% end

