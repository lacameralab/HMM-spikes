function [col2] = plotviterbi(firings,stateseq,N,timev,varargin)
% [col2] = plotviterbi(firings,pstates,N,timev,varargin)
% This function plots the state sequence decoded by Viterbi algorith.
% This plot will shift towards left side because it plots the pstates at
% the left boundary of the time window (dt, which was used to generate
% spike counts for training.)
% 
% Inputs:
%   firings: two column matrix storing all the recorded neurons' spike time
%            in trial k. The first column is the spike time and the second
%            column is the index of the neuron which fires at the
%            corresponding time on the same row. 
%   stateseq: vector, row i is the probability of being in state i at the
%            corresponding time bin.
%   N: number of neurons
%   timev: time vector
%   'isnewfig': if isnewfig=0, plot on the existed figure; else, generate
%               new figure. Default: 1.
%   'fontsize': font size of x, y label. Defautl: 18.
%   'spikeshape': string, the shape of each spike, 'dot' or 'vertline'
%   'statecolor': matrix, row i is the color vector for state i
% 
% Outputs:
%   'col2': color matrix used in plot
% 
% Tianshu Li

[varargin{:}] = convertStringsToChars(varargin{:});

% isnewfig
isisnewfig = cellfun(@(x) strcmp(x,'isnewfig'), varargin);
idxisnewfig = find(isisnewfig==1, 1, 'first');
if idxisnewfig
    isnewfig = varargin{idxisnewfig+1};
end

% fontsize
isfontsize = cellfun(@(x) strcmp(x,'fontsize'), varargin);
idxfontsize = find(isfontsize==1, 1, 'first');
if idxfontsize
    fntsz = varargin{idxfontsize+1};
end

% spikeshape, the shape of spikes, dot or vertline
isspikeshape = cellfun(@(x) strcmp(x,'spikeshape'), varargin);
idxspikeshape = find(isspikeshape==1, 1, 'first');
if idxspikeshape
    spikeshape = varargin{idxspikeshape+1};
else
    spikeshape = 'dot';
end

% statecolor, the color matrix for each state
isstatecolor = cellfun(@(x) strcmp(x,'statecolor'), varargin);
idxstatecolor = find(isstatecolor==1, 1, 'first');
if idxstatecolor
    setcol = varargin{idxstatecolor+1};
end


if ~exist('isnewfig','var')
    isnewfig = 1;
end

%  -------------   Set color   ---------------
m = size(unique(stateseq));
if ~exist('setcol','var')
    cols = lines;
    setcol = cols;
    deltacol = 7; % lines is a vector with 7 colors in loop
    for i = 1:m
        setcol(i*deltacol+1:(i+1)*deltacol,:) = mod(cols(i*deltacol+1:(i+1)*deltacol,:)+rand(deltacol,3),1);
    end
end
% --------------------------------------------

starttime = timev(1);
endtime = timev(end);
dt = timev(2)-timev(1);
timev_ori = timev;
timev = (timev-starttime)/dt;
timev = round(timev);

stateseq = [stateseq, stateseq(:,end)];

if isnewfig
    figure; clf;
end

if ~exist('fntsz','var')
    fntsz=18;
end

hold all; 

nState = max(stateseq);

% use color to identify states
transtimeind = [1, find(diff(stateseq))+1, length(timev)];
transtimeind = sort(transtimeind);
for i=1:length(transtimeind)-1
    x1 = [transtimeind(i),transtimeind(i+1)];
    x2 = [x1, fliplr(x1)];
    s = stateseq(transtimeind(i)); % current state
    if ~isempty(x2)
        patch(timev_ori(x2), [0,0,N+1,N+1], setcol(s,:),'FaceAlpha',0.5,'EdgeColor','none');
    end
end

col2 = setcol(1:m,:);

rasterPlot(firings,'newfigure',0,'spikeshape',spikeshape,'markersize',2); % raster at top (on top of state color patch)

hold off;
figset(gca,'time (sec)','neuron index',fntsz);

end