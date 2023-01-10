function [] = rasterPlot(spkTrains, varargin)
% =======================================================================
% [] = rasterPlot(spkTrains, varargin): 
%      generate a raster plot of input spike trains. 
%
% Parameters:
%   spkTrains: 
%     Data type 1: cell
%     Description: the nth element of the cell is a vector storing spike
%     time of the nth neuron
%     Data type 2: plot type
%     Description: a 2-column matrix, the first column is the firing time
%     of all the neurons, the second column is the corresponding neuron
%     index of the neuron which fires at that time (the firing time on the
%     same row, 1st column in the matrix) (or the nth trial) . We use this
%     data type to plot raster plot.
%   newfigure:
%     If newfigure, create a new figure, else, plot on the current figure.
%     Default value is 1.
%   color:
%     Data type: vector with 3 or 4 entries
%     Description: color vector of rasters, default: [0,0,0] (black)
%   stateSeq:
%     Data type: matrix
%     Description: a 2 row matrix, the first row is the time when state
%     changes, the second row is the state transit to. 
%   fontsize:
%     Datatype: number
%     Description: font size
%   linewidth:
%     Datatype: number
%     Description: dot size
%   'spikeshape':
%     Datatype: string
%     Description: the shape of each spike, 'dot' or 'vertline'
%   'markersize':
%     Datatype: number
%     Description: the size of the spike marker. Marker size of dot or
%     linewidth of raster (short vertical bar).
% =======================================================================

% If stateSeq is not given, do not visualize states.
if ~exist('newfigure','var')
    newfigure = 1;
end

is_state = 0;

if ~exist('col','var')
    col = [0,0,0];
end

[varargin{:}] = convertStringsToChars(varargin{:});

% state sequence
isstateSeq = cellfun(@(x) strcmp(x,'stateSeq'), varargin);
idxstateSeq = find(isstateSeq==1, 1, 'first');
if idxstateSeq
    stateSeq = varargin{idxstateSeq+1};
    is_state = 1;
end

% line width
islinewidth = cellfun(@(x) strcmp(x,'linewidth'), varargin);
idxlinewidth = find(islinewidth==1, 1, 'first');
if idxlinewidth
    linewidth = varargin{idxlinewidth+1};
end

% color
iscolor = cellfun(@(x) strcmp(x,'color'), varargin);
idxcolor = find(iscolor==1, 1, 'first');
if idxcolor
    col = varargin{idxcolor+1};
end

% generate new figure
isnewfigure = cellfun(@(x) strcmp(x,'newfigure'), varargin);
idxnewfigure = find(isnewfigure==1, 1, 'first');
if idxnewfigure
    newfigure = varargin{idxnewfigure+1};
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

% marker size
ismarkersize = cellfun(@(x) strcmp(x,'markersize'), varargin);
idxmarkersize = find(ismarkersize==1, 1, 'first');
if idxmarkersize
    markersize = varargin{idxmarkersize+1};
end

if ~exist('markersize','var')
    if strcmp(spikeshape,'dot')
        markersize=18;
    elseif strcmp(spikeshape,'vertline')
        markersize=1;
    end
end


if ~exist('linewidth','var')
    linewidth=18;
end

if ~exist('fntsz','var')
    fntsz=18;
end

if isempty(spkTrains)
    if newfigure
        figure; hold on;
        set(gcf, 'position', [254   175   894   508]);
    else
        hold on;
    end
    return;
end

% Check the data type of spkTrains, and change the data type to plot type.
if iscell(spkTrains)
    if length(spkTrains)==1 && isempty(spkTrains{:})
        if newfigure
            figure; hold on;
            set(gcf, 'position', [254   175   894   508]);
        else
            hold on;
        end
        return;
    end
    spkTrains = spkTrains2plottype(spkTrains);
end


% Plot
if newfigure
    figure; hold on;
    set(gcf, 'position', [254   175   894   508]);
else
    hold on;
end
if isempty(spkTrains)
    return;
end
if strcmp(spikeshape,'dot')
    plot(spkTrains(:,1),spkTrains(:,2),'.','color',col,'linewidth',markersize);
elseif strcmp(spikeshape,'vertline')
    plot([spkTrains(:,1),spkTrains(:,1),NaN(size(spkTrains(:,1)))]',[spkTrains(:,2)-0.2,spkTrains(:,2)+0.2,NaN(size(spkTrains(:,2)))]','-','color',col,'linewidth',markersize);
end
xlim([min(spkTrains(:,1))-0.1 max(spkTrains(:,1))+0.1]);
N = max(spkTrains(:,2));
ylim([0 N+1]);
% only label the neuron index in yaxis (remove the N+1 yticklabel)
yticks_rec = yticks;
yticks_rec = round(yticks_rec);
yticks_rec = unique(yticks_rec);
if yticks_rec(end) > N
    yticks_rec = yticks_rec(1:end-1);
end
yticks(yticks_rec);
% yticks(1:N);
fntsize = 20;
title('Raster plot');
figset(gca,'time (sec)','Trial index',fntsize);


% If state sequence (stateSeq) is given, visualize the state transition.
if is_state
    for i = 1:size(stateSeq,2)
        plot([stateSeq(1,i) stateSeq(1,i)], ylim, 'k');
        text(stateSeq(1,i),max(ylim),strcat('s',string(stateSeq(2,i))),'color','r','fontsize',20);
    end
    ylabel('Neuron index');
end
if newfigure
    hold off;
end

end