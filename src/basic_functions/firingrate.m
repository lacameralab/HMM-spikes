function [y_smoothed,x] = firingrate(firings,neuind,dt,edges)
% [y_smoothed,x] = firing_rate(firings,neuind,dt,edges)
% This function uses Gaussian kernel to smooth the spike train (bandwidth = dt (sec)).
% Inputs:
%   firings: a 2-column matrix, the first column is the firing time
%     of all the neurons, the second column is the corresponding neuron
%     index of the neuron which fires at that time (the firing time on the
%     same row, 1st column in the matrix) (or the nth trial).
%   neuind: interger, the index of neuron of interest
%   dt: time bin, same unit as firings
%   edges: bin edges
% Outputs:
%   y_smoothed: smoothed firing rate matrix, each row is a neuron
%   x: time vector for y_smoothed
% 
% Tianshu Li

if ~exist('dt','var')
    if exist('edges','var')
        dt = edges(2)-edges(1);
    else
        dt = 1;
    end
end

y = firings(firings(:,2)==neuind,1);

if isempty(y)
    if exist('edges','var')
        y_smoothed = zeros(1,length(edges));
    else
        y_smoothed = 0;
    end
else
    y_estimator = fitdist(y,'Kernel','BandWidth',dt);
    if exist('edges','var')
        x = edges;
    else
        x = min(firings(:,1)):0.1:max(firings(:,1));
    end
    y_smoothed = pdf(y_estimator,x)*length(y);
end

end