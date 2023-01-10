function [seq,newfirings] = firings2seq(firings,neuind,timev)
% [seq,newfirings] = firings2seq(firings,neuind,timev)
%       This function is used to transfer "firings" to obsevation sequence.
% Inputs:
%   firings: a two column matrix storing all the recorded
%      neurons' spike time in trial k. The first column is the spike time
%      and the second column is the index of the neuron which fires at the
%      corresponding time on the same row.
%   neuind: chosen neurons' indices
%   timev: time vector, the intervals are the time windows for each
%      observation.
% Outputs:
%   seq: matrix, observation sequences, each row is one trial
%   newfirings: the spike data after removing the duplicated spikes in each
%      time bin. This one is used to generate the observation sequences.
% 
% Notice: firings is in second, timev and dt's units are second too.

dt = timev(2)-timev(1);
timev = 1000*timev; % ms
dt = 1000*dt; % ms
firings(:,1) = firings(:,1)*1000; % sec to ms
% at most one spike in each time bin
[newfr,newfirings] = preparedata4HMM(firings, timev, neuind, dt, 1);
seq = 1;
for k = 1:length(neuind)
    i = neuind(k);
    seq = seq + i*(newfr(i,:)==1);
end

end