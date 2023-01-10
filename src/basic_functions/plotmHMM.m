function [logpseq, col2] = plotmHMM(seq, estTRANS, estEMIS, firings, timev, varargin)
% [logpseq] = plotmHMM(seq, estTRANS, estEMIS, firings, timev, varargin)
% This function is used to plot the decoding of mHMM on ONE trial. It will
% reture the loglikelihood of the HMM model on this trial - "logpseq".
% Inputs:
%   seq: symbol sequence, the input to mHMM (one row)
%   estTRANS: estimated transition matrix
%   estEMIS: estimated emission probability matrix
%   firings: the spiking data of the trial, 2 columns, first column is the
%       spike time, second column is the neuron which spikes at the
%       corresponding time.
%   timev: time vector, corresponds to "seq"
%   'stateSeq': state sequence, two rows, first row is the transition time,
%       second row is the state transits to.
%   'isnewfig': if 1, open a new figure, otherwise, plot on the existed
%       figure, default: 1
%   'fontsize': integer, font size, default: 18
%   'statecolor': matrix, row i is the color vector for state i
%   'spikeshape': string, the shape of each spike, 'dot' or 'vertline'
% Outputs:
%   'logpseq': loglikelihood sequence across time
%   'col2': color matrix used in plot
% 
% Tianshu Li
% Dec. 2020

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

if ~exist('fntsz','var')
    fntsz=18;
end

%  -------------   Set color   ---------------
m = size(estTRANS,1);
if ~exist('setcol','var')
    cols = lines;
    setcol = cols;
    deltacol = 7; % lines is a vector with 7 colors in loop
    for i = 1:m
        setcol(i*deltacol+1:(i+1)*deltacol,:) = mod(cols(i*deltacol+1:(i+1)*deltacol,:)+rand(deltacol,3),1);
    end
end
% --------------------------------------------

%% visualization
timev_ori = timev;
dt = timev(2)-timev(1);
starttime = timev(1);
endtime = timev(end);
timev = (timev-starttime)/dt;
timev = round(timev);
if isnewfig
    figure; clf;
end
hold on;

% Probability of each state
[pstates,logpseq] = hmmdecode(seq(1,:),estTRANS,estEMIS);

nState = size(estTRANS,1);
N = max(unique(firings(:,2)));

col2 = [];
yyaxis right;
for i = 1:nState
    h = plot(starttime+timev(1:end-1)*dt,pstates(i,:),'-','color',setcol(i,:),'linewidth',3);
    col2 = [col2; get(h,'Color')];
end
figset(gca,'time (sec)','Probability in each state',fntsz);

yyaxis left;
% Take the state where its pstate is not smaller than 80% as the state at
% that time. If there is no state meet this condition, do not assign any
% state at that time.
% plot the decoded states (those with max prob in each bin): 
[maxP, indMax]=max(pstates);
ind=find(maxP>0.8); % keep only states inferred with high prob.
% use color to identify states
transtimeind = [1, find(indMax(2:end)-indMax(1:end-1)~=0), length(timev)];
for i=2:length(transtimeind)
    x = timev(transtimeind(i-1):transtimeind(i));
    x = x(ismember(x,ind));
    xjump = [0,find(x(2:end)-x(1:end-1)~=1),length(x)];
    for j = 1:length(xjump)-1
        x1 = x(xjump(j)+1:xjump(j+1));
        x2 = [x1, fliplr(x1)];
        s = indMax(transtimeind(i-1)+1); % current state
        inBetween = [(N+1)*[pstates(s,x1)>0.8], fliplr(0*x1)]; % colored areas are rectangles covering the time bins in that state
        if ~isempty(x2)
            patch(timev_ori(x2), inBetween', col2(s,:),'FaceAlpha',0.5,'EdgeColor','none');
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

end