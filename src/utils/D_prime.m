function D_prime = D_prime(X1, X2)
% This function computes the D' metric between the 2 distributions that X1
% and X2 follow.
% The D' (D-prime) metric is a measure used in signal detection theory
% (SDT) to quantify the difference between two distributions, often
% representing signal and noise.
% 
% Inputs:
%   X1, X2: data samples of two variables.
% 
% Outputs:
%   D_prime: we take the absolute value of D' in this function.
% 
% Tianshu Li
% Mar. 31, 2025

n1 = length(X1);
n2 = length(X2);

mu1 = mean(X1);
mu2 = mean(X2);

sigma1 = std(X1);
sigma2 = std(X2);

sigma_pooled = sqrt(((n1-1)*sigma1^2+(n2-1)*sigma2^2)/(n1+n2-2)); % pooled standard deviation
D_prime = abs(mu1-mu2)/sigma_pooled;

end