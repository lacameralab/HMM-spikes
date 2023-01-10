function [logalpha, logbeta] = phmm_lalphabeta(m,x,lambda,Gamma,delta,lP_xt)
% [logalpha, logbeta] = phmm_lalphabeta(m,x,lambda,Gamma,delta,lP_xt)
% Calculate the log of forward (alpha) and backward (beta) probabilities
% Inputs: 
%   m: number of states
%   x: observation matrix
%   lambda: firing rate matrix
%   Gamma: transition matrix
%   delta: stationary distribution (optional)
%   lP_xt: log(P(xt)), P(xt): probability of x=xt at time t, p(i,t) = P(xt|Ci) (optional)
% Outputs:
%   logalpha: log(alpha), alpha is the forward probability
%   logbeta: log(beta), beta is the backward probability
% 
% Tianshu Li - edited from the R code in Zucchini's book A.2.3.

if nargin < 5
    delta = ones(1,m)/(eye(m)-Gamma+ones(m));
end

n = size(x,2); % length of observations (time steps)

% logalpha
% ------------------------------
logalpha = zeros(m,n);
if nargin < 6
    % log(P(xt))
    lP_xt = zeros(m,n); % P(xt), p(i,t) = P(xt|Ci)
    for i = 1:m
        lP_xt(i,:) = sum(-lambda(:,i)+log(lambda(:,i)).*x-log(factorial(x)),1);
    end
end
P_xt = exp(lP_xt);
alpha = delta.*P_xt(:,1)';
logscale = log(sum(alpha));
alpha = alpha/sum(alpha);
logalpha(:,1) = (log(alpha)+logscale)';
for t = 2:n
    alpha = (alpha*Gamma).*P_xt(:,t)';
    logscale = logscale + log(sum(alpha));
    alpha = alpha/sum(alpha);
    logalpha(:,t) = (log(alpha)+logscale)';
end
% -------------------------------
logbeta = zeros(m,n);

beta = ones(1,m)/m;
logscale = log(m);
for t = n-1:-1:1
    beta = Gamma*(P_xt(:,t+1).*beta');
    beta = beta';
    logbeta(:,t) = (log(beta)+logscale)';
    logscale = logscale + log(sum(beta));
    beta = beta/sum(beta);
end

end