function [spkc,neuind] = spikecount(firings, dt, timev, neuind)
% [spkc] = spikecount(firings, dt, timev, neuind)
% Count spikes in each time bin.
% Inputs:
%   firings: 
%       Data type 1: plotType matrix (2 columns, first column is the firing time,
%            second column is the corresponding neuron's index that fires at that
%            time.)
%       Data type 2: cell, the nth element of the cell is a vector storing
%            spike time of the nth neuron .
%   dt: time window for calculating firing rate, default: [] (count all spikes in the trial)
%   timev: time vector for window edges, default: [-Inf, Inf] (count all spikes in the trial)
%   neuind: index of neurons, default is the all the unique numbers in the
%           second column of firings (sorted).
% Return: 
%    spkc: a nNeu*nTimewindow matrix, row i is the firing rate in each
%        time window of neuron i.

% (nNeu: number of neurons, usually, nNeu=N, but when nNeu=Ne, this function
%       only caculate E neurons' firing rates. not used anymore)

if nargin < 3
    if nargin < 2
        dt = [];
        timev = [-Inf, Inf];
    else
        T = max(firings(:,1));
        timev = 0:dt:T;
        if timev(end) ~= T
            timev = [timev, T];
        end
    end
end

% neuind = unique(firings(:,2)); % neuron index that fires
% if nargin < 4
%     nNeu = max(neuind); % number of neurons
% end
% neuind(neuind>nNeu) = []; % remove the neurons whose index is larger than nNeu (this is used when only wants to caculate E neurons' firing rates)

% Check the data type of spkTrains, and change the data type to plot type.
if iscell(firings)
    if length(firings)==1 && isempty(firings{:})
        firings = [];
    else
        firings = spkTrains2plottype(firings);
    end
end
if isempty(firings)
    spkc = [];
    return;
end
if nargin < 4
    neuind = unique(firings(:,2)); % index of neurons that fires
    neuind = sort(neuind);
end


% flag = 0; % if flag = 1, the last time window is shorter than dt
% timev = 0:dt:T;
% if timev(end) ~= T
%     timev = [timev, T];
% %     flag = 1;
% end

nNeu = length(neuind);
spkc = zeros(nNeu, length(timev)-1);

for i = 1:nNeu
    spkc(i,:) = histcounts(firings(firings(:,2)==neuind(i),1),timev);
%     ind = find(firings(:,2)==i);
%     for j = 1:length(timev)-2
%         fr(i,j) = length(find(firings(ind,1)<timev(j+1)&firings(ind,1)>=timev(j)))/dt;
%     end
%     j = j + 1;
%     if flag == 1
%         fr(i,j) = length(find(firings(ind,1)<timev(j+1)&firings(ind,1)>=timev(j)))/(T-timev(end-1));
%     else
%         fr(i,j) = length(find(firings(ind,1)<timev(j+1)&firings(ind,1)>=timev(j)))/dt;
%     end
end

end