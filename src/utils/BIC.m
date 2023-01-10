function [bic] = BIC(m,N,Nbin,ntrial,LL,HMMtype)
% [bic] = BIC(m,N,Nbin,ntrial,LL,HMMtype)
% Inputs:
%   m: number of states
%   N: number of neurons
%   Nbin: number of observations (time bin) for one trial, put T/dt if T is
%       the total time of one trial, dt is bin size. 
%   ntrial: number of trials used to calculated LL
%   LL: log-likelihood of training set
%   HMMtype: type of HMM, 'bHMM' or 'pHMM' for binomial HMM and Poisson
%       HMM.
% Ouput:
%   bic: BIC of the model

LL = real(LL);
if strcmp(HMMtype,'pHMM')
    bic = -2*LL+(m*(m-1)+m*N)*log(ntrial*N*Nbin);
elseif strcmp(HMMtype,'bHMM')
    bic = -2*LL+(m*(m-1)+m*N)*log(ntrial*Nbin);
end
end
