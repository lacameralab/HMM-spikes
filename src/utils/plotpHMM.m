function [col2] = plotpHMM(firings,pstates,N,timev,varargin)
% [col2] = plotpHMM(firings,pstates,N,timev,varargin)
% This plot will shift towards left side because it plots the pstates at
% the left boundary of the time window (dt, which was used to generate
% spike counts for training.)
% 
% Inputs:
%   firings: two column matrix storing all the recorded neurons' spike time
%            in trial k. The first column is the spike time and the second
%            column is the index of the neuron which fires at the
%            corresponding time on the same row. 
%   pstates: matrix, row i is the probability of being in state i at the
%            corresponding time bin.
%   N: number of neurons
%   timev: time vector
%   'stateSeq': state sequence
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

% state sequence
isstateSeq = cellfun(@(x) strcmp(x,'stateSeq'), varargin);
idxstateSeq = find(isstateSeq==1, 1, 'first');
if idxstateSeq
    stateSeq = varargin{idxstateSeq+1};
end

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
m = size(pstates,1);
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

pstates = [pstates, pstates(:,end)];

if isnewfig
    figure; clf;
end

if ~exist('fntsz','var')
    fntsz=18;
end

hold all; 

nState = size(pstates);

yyaxis left;
[maxP, indMax]=max(pstates);
ind=find(maxP>0.8); % keep only states inferred with high prob.
% use color to identify states
transtimeind = [0, find(diff(indMax))+1, ind(diff(ind)>1)+1, length(timev)];
transtimeind = sort(transtimeind);
for i=2:length(transtimeind)
    x = timev(transtimeind(i-1)+1:transtimeind(i));
    x = x(ismember(x,ind));
    xjump = [0,find(diff(x)~=1)+1,length(x)];
    for j = 1:length(xjump)-1
        x1 = x(xjump(j)+1:xjump(j+1));
        x2 = [x1, fliplr(x1)];
        s = indMax(transtimeind(i)-1); % current state
        inBetween = [(N+1)*[pstates(s,x1)>0.8], fliplr(0*x1)]; % colored areas are rectangles covering the time bins in that state
        if ~isempty(x2)
            patch(timev_ori(x2), inBetween', setcol(s,:),'FaceAlpha',0.5,'EdgeColor','none');
        end
    end
end

yyaxis left;
rasterPlot(firings,'newfigure',0,'spikeshape',spikeshape,'markersize',2); % raster at top (on top of state color patch)

if exist('stateSeq','var')
    stairs([stateSeq(1,:),endtime],[stateSeq(2,:),stateSeq(2,end)]-0.8,'k','linewidth',3);
    ylim([0,N+1]);
end
hold off;
figset(gca,'time (sec)','neuron index',fntsz);
title(sprintf('PHMM (dt = %g ms)',dt*1000));

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

% plot the state probabilities
col2 = [];
yyaxis right; hold on;
ylim([0,1]);
for i = 1:nState
    h = plot(starttime+timev*dt,pstates(i,:),'-','color',setcol(i,:),'linewidth',3);
    col2 = [col2; get(h,'Color')];
end

figset(gca,'time (sec)','Probability in each state',fntsz);

end
