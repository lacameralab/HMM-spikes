function [stateprobs] = phmm_stateProbs(x,m,lambda,Gamma,delta)
% [stateprobs] = phmm_stateProbs(x,m,lambda,Gamma,delta)
% This function calculates the probability of each state at each time bin.
% Inputs: 
%   x: observation matrix
%   m: number of states
%   lambda: firing rate matrix
%   Gamma: transition probability matrix
%   delta: initial state distribution (optional)
% Outputs:
%   stateprobs: matrix, each row is the probability of staying in one state at a
%       time bin.
% 
% - Tianshu Li

if nargin < 5
    delta = abs(ones(1,m)/(eye(m)-Gamma+ones(m))); % set initial state distribution to stationary distribution
    % delta = zeros(1,m); delta(1) = 1; % assume every trial starts at state 1
end

n = size(x,2); % number of time bins
[la, lb] = phmm_lalphabeta(m,x,lambda,Gamma,delta); % forward and backward probabilities: log(alpha), log(beta)

c = max(la(:,n));
llk = c+log(sum(exp(la(:,n)-c))); % loglikelihood
stateprobs = zeros(m,n);
for i = 1:n
    stateprobs(:,i) = exp(la(:,i)+lb(:,i)-llk);
end

end
