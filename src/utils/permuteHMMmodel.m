function [permuted_HMM] = permuteHMMmodel(HMMmodel)
% shuffle/permute the firing rate of each neuron across states to generate
% samples.
% Input:
%   HMMmodel: structure, HMMmodel.Gamma is the transition probability
%   matrix and HMMmodel.Lambda is the firing rate matrix
% 
% Output:
%   permuted_HMM: structure, same as HMMmodel, the permuted HMM model
% 
% Tianshu Li
% Apr. 10, 2025

[N, m] = size(HMMmodel.Lambda);
permuted_HMM.Lambda = 0*HMMmodel.Lambda;
permuted_HMM.Gamma = HMMmodel.Gamma;
% shuffle/permute the firing rate of each neuron across states to generate
% samples (do it on bestmodel)
for i = 1:N
    permuted_HMM.Lambda(i,:) = HMMmodel.Lambda(i, randperm(m));
end
% permute the off diagonal values in Gamma
for i = 1:m
    permuted_HMM.Gamma(i,setdiff(1:m,[i])) = HMMmodel.Gamma(i, randperm(m-1));
end