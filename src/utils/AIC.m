function [aic] = AIC(m,N,LL)
% [aic] = AIC(m,N,LL)
% 
% Inputs:
%   m: number of states
%   N: number of neurons
%   LL: log-likelihood of training set
% Ouput:
%   aic: AIC of the model

aic = -2*LL+2*(m*(m-1)+m*N);

end