function [newfr,newfirings] = preparedata4HMM(firings, timev, neuind, dt, isnewfirings)
% [newfr,newfirings] = preparedata4HMM(firings, timev, neuind, dt, isnewfirings)
% This function is used for multinoulli HMM (mHMM).
% Only deal with selected neurons (neuind) in HMM.
% Match the HMM model to the data segmented in 1 ms bins. Given the
% short duration of the bins, we assumed that in each bin either zero or
% one spikes are emitted by each neuron. Neglect the simultaneous firing of
% two or more neurons (a rare event): if more than one neuron fired in a
% given bin, a single spike was randomly assigned to one of the firing
% neurons. This function is used to prepare data for HMM to meet the
% previous assumption.
% Parameters:
%   firings - first column is the firing time, second column is the
%             corresponding neuron that fires.
%   timev - time vector.
%   neuind - vector, selected neurons' indices
%   dt - time bin for HMM (same unit as firings), default is 1.
%   isnewfirings - bool or interger, if isnewfirings is true (or non zero
%                  number), calculate newfirings.
% Output:
%   newfr - N*nBin matrix, number of spikes (zero or one) in each time bin
%           for each neuron (if more than one, randomly take one of them).
%   newfirings - 2 column matrix, the firings matrix after removing spikes
%                where there are more than one spikes in that time bin.

if nargin < 5
    isnewfirings = 1;
end
if nargin < 4
    dt = 1; % same unit as firings
end

% Step 1: remove neurons not selected.
newfirings = firings;
indE = find(~ismember(firings(:,2),neuind));
newfirings(indE,:) = [];

% Step 2: randomly choose one spikes in each bin when there are more than one
% spikes and remove others for each neuron.
fr = spikecount(firings,dt,timev,neuind);
newfr = fr;
indmorespkfr = find(fr>1); % index of each neuron that spike more than one time in one time bin
newfr(indmorespkfr) = 1;

if length(neuind) > 1
    % Step 3: randomly choose one neuron in each bin fires when there are more
    % than one neuron spike in that bin.
    indmulneuspk = find(sum(newfr)>1); % index of time bin where there are more than one neuron spikes
    for i=indmulneuspk
        indneu = find(newfr(:,i)>0); % neuron indexs that fires in time bin i
        chosenneuind = indneu(randi(length(indneu))); % randomly choose one neuron
        newfr(:,i)=0; % assign all neurons in time bin i have 0 spike
        newfr(chosenneuind,i)=1; % assign the chosen neuron 1 spike in time bin i
    end
end

if isnewfirings
    % reset data s.t. each bin only has zero or one spike.
    for i = 1:length(neuind)
        indmorespk = find(fr(i,:)>1);
        for j = indmorespk
            ind = find(newfirings(:,2)==neuind(i)&newfirings(:,1)<timev(j+1)&newfirings(:,1)>=timev(j));
            randind = randi(length(ind)); % randomly choose one spike and remove the others
            ind(randind) = [];
            newfirings(ind,:) = [];
        end
    end
end
end