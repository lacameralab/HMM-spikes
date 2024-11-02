function [pStates,logpSeq] = phmm_decode(x,Gamma,lambda)
% [pStates,logpSeq] = phmm_decode(x,Gamma,lambda,varargin)
% 
% phmm_decode calculates the posterior state probabilities of a spike count
% sequence given a Poisson HMM model defined by transition probability
% matrix Gamma and firing rate matrix lambda.
% N is the number of neurons, T is the number of time bins, and m is the
% number of states.
% Inputs:
%   x: spike count sequence (NxT)
%   Gamma: transition probability matrix (mxm)
%   lambda: firing rate matrix (Nxm), the unit is number of spikes in a
%       time bin. 
% Outputs:
%   pStates: posterior state probabilities
%   logpSeq: log of the probability of the whole sequence
% 
% Modified code from Matlab's function HMMDECODE.

numStates = size(Gamma,1);
checkGamma = size(Gamma,2);
if checkGamma ~= numStates
    error('Transition probability matrix (Gamma) must be square.');
end

% number of columns of lambda must be same as number of states
checklambda  = size(lambda,2);
if checklambda ~= numStates
    error('The firing rate matrix (lambda) must have the same number of states as the transition probability matrix (Gamma).');
end

% add extra data to start to make algorithm cleaner at f0 and b0
x = [zeros(size(lambda,1),1), x ];
T = size(x,2);

% P(xt)
P_xt = zeros(numStates,T); % P(xt), p(i,t) = P(xt|Ci)
for i = 1:numStates
    P_xt(i,:) = exp(sum(-lambda(:,i)+log(lambda(:,i)).*x-log(factorial(x)),1));
end

% introduce a scaling factor for alpha (forward)
fs = zeros(numStates,T); % scaled alpha
fs(1,1) = 1;  % assume that we start in state 1.
s = zeros(1,T);
s(1) = 1;
for count = 2:T
    for state = 1:numStates
        fs(state,count) = P_xt(state,count) .* (sum(fs(:,count-1) .*Gamma(:,state)));
    end
    % scale factor normalizes sum(fs,count) to be 1. 
    s(count) =  sum(fs(:,count));
    fs(:,count) =  fs(:,count)./s(count);
end

% use the scale factor for beta (backward)
bs = ones(numStates,T); % scaled beta
for count = T-1:-1:1
    for state = 1:numStates
      bs(state,count) = (1/s(count+1)) * sum( Gamma(state,:)'.* bs(:,count+1) .* P_xt(:,count+1)); 
    end
end

%  The  actual backward and  probabilities can be recovered by using
%  scales = cumprod(s, 'reverse'); 
%  b = bs.*repmat([scales(2:end), 1],size(bs,1),1);

logpSeq = sum(log(s));
pStates = fs.*bs;

% get rid of the column that we stuck in to deal with the f0 and b0 
pStates(:,1) = [];

end

