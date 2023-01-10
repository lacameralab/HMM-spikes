function [stateprobs] = phmm_stateProbs(x,m,lambda,Gamma,delta)
% [stateprobs] = phmm_stateProbs(x,m,lambda,Gamma,delta)
% This function calculates the probability of each state at each time bin.
% Inputs: 
%   x: observation matrix
%   m: number of states
%   lambda: firing rate matrix
%   Gamma: transition matrix
%   delta: stationary distribution (optional)
% Outputs:
%   stateprobs: matrix, each row is the probability of one state at each
%       time bin.
% 
% Tianshu Li

if nargin < 5
    delta = abs(ones(1,m)/(eye(m)-Gamma+ones(m)));
end

n = size(x,2);
[la, lb] = phmm_lalphabeta(m,x,lambda,Gamma,delta);

c = max(la(:,n));
llk = c+log(sum(exp(la(:,n)-c)));
stateprobs = zeros(m,n);
for i = 1:n
    stateprobs(:,i) = exp(la(:,i)+lb(:,i)-llk);
end

end