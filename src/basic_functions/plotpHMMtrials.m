function [ah] = plotpHMMtrials(spkc,lambda,Gamma,timev,varargin)
% [ah] = plotpHMMtrials(spkc,lambda,Gamma,timev,varargin)
% This function plots the decoding of all trials in the same plot.
% Inputs:
%   spkc: spike count cell array, each cell contains the spikc count matrix
%       of one trial.
%   lambda: firing rate matrix, unit: inverse of the unit of timev
%   Gamma: transition probability matrix
%   timev: time vector
%   varargin:
%       'isnewfig': if 1, generate a new figure to plot on
%       'isfiglabels': if 1, add labels to the plot. Set this variable to 0
%           for presentation or paper. Default: 1
%       'fontsize': font size, default: 18
%       'isplotfr': if 1, plot firing rate matrix, default: 1
%       'statecolor': matrix, row i is the color vector for state i
%       'stateshape': 'square' or 'probcurve' 
%           If 'square', each state is colored as a square between the trial
%           indices. 
%           If 'probcurve', only the area under the probability curve
%           of the estimated state will be colored.
%           Default: 'square'
% 
% Tianshu Li
% Jan 27, 2021

if ~exist('phmm_stateProbs') addpath '../pHMM'; end

nTrial = length(spkc);
% --------------  Parameters  ----------------
[varargin{:}] = convertStringsToChars(varargin{:});

% isnewfig
isisnewfig = cellfun(@(x) strcmp(x,'isnewfig'), varargin);
idxisnewfig = find(isisnewfig==1, 1, 'first');
if idxisnewfig
    isnewfig = varargin{idxisnewfig+1};
else
    isnewfig = 1;
end

% isfiglabels
isisfiglabels = cellfun(@(x) strcmp(x,'isfiglabels'), varargin);
idxisfiglabels = find(isisfiglabels==1, 1, 'first');
if idxisfiglabels
    isfiglabels = varargin{idxisfiglabels+1};
else
    isfiglabels = 1;
end

% fontsize
isfontsize = cellfun(@(x) strcmp(x,'fontsize'), varargin);
idxfontsize = find(isfontsize==1, 1, 'first');
if idxfontsize
    fntsz = varargin{idxfontsize+1};
else
    fntsz=18;
end

% isplotfr, if 1, plot firing rate matrix
isisplotfr = cellfun(@(x) strcmp(x,'isplotfr'), varargin);
idxisplotfr = find(isisplotfr==1, 1, 'first');
if idxisplotfr
    isplotfr = varargin{idxisplotfr+1};
else
    isplotfr = 1;
end

% statecolor, the color matrix for each state
isstatecolor = cellfun(@(x) strcmp(x,'statecolor'), varargin);
idxstatecolor = find(isstatecolor==1, 1, 'first');
if idxstatecolor
    setcol = varargin{idxstatecolor+1};
end

% stateshape, the shape of the presentation of states
isstateshape = cellfun(@(x) strcmp(x,'stateshape'), varargin);
idxstateshape = find(isstateshape==1, 1, 'first');
if idxstateshape
    stateshape = varargin{idxstateshape+1};
else
    stateshape = 'square';
end

% --------------   Set color   ---------------
m = size(Gamma,1);
N = size(lambda,1);
if ~exist('setcol','var')
    cols = lines;
    setcol = cols;
    deltacol = 7; % lines is a vector with 7 colors in loop
    for i = 1:m
        setcol(i*deltacol+1:(i+1)*deltacol,:) = mod(cols(i*deltacol+1:(i+1)*deltacol,:)+rand(deltacol,3),1);
    end
end

% --------------- Preparation ----------------
% when all the trials have the same timev
starttime = timev(1);
endtime = timev(end);
dt = timev(2)-timev(1);
timev_ori = timev;
timev = (timev-starttime)/dt;
timev = round(timev);

% ------------------  Plot  ------------------
if isnewfig
    figure; clf;
    set(gcf,'position',[193         404        1292         420]);
end
ah = subplot(3,m,[1:2*m]);
hold all; 
for k = 1:nTrial
    pstates = phmm_stateProbs(spkc{k},m,lambda*dt,Gamma);
    pstates = [pstates, pstates(:,end)];
    
    [nState, ntimebin] = size(pstates);
    
    [maxP, indMax]=max(pstates);
    ind=find(maxP>0.8); % keep only states inferred with high prob.
    
    transtimeind = [0, find(diff(indMax))+1, ind(diff(ind)>1)+1, length(timev)];
    transtimeind = sort(transtimeind);
    for i=2:length(transtimeind)
        x = timev(transtimeind(i-1)+1:transtimeind(i));
        x = x(ismember(x,ind));
        xjump = [0,find(diff(x)~=1)+1,length(x)];
        for j = 1:length(xjump)-1
            x1 = x(xjump(j)+1:xjump(j+1));
            if ~isempty(x1) && x1(end) < ntimebin
                x1 = [x1, x1(end)+1];
            end
            x2 = [x1, fliplr(x1)];
            s = indMax(transtimeind(i)-1); % current state
            inBetween = [k+pstates(s,x1), k+fliplr(0*x1)];
            if ~isempty(x2)
                if strcmp(stateshape,'square')
                    xblock = [min(timev_ori(x2)),max(timev_ori(x2)),max(timev_ori(x2)),min(timev_ori(x2))];
                    patch(xblock, [k,k,k+1,k+1], setcol(s,:),'FaceAlpha',0.5,'EdgeColor','none');
                elseif strcmp(stateshape,'probcurve')
                    patch(timev_ori(x2), inBetween', setcol(s,:),'FaceAlpha',0.5,'EdgeColor','none');
                end
            end
        end
    end
end

xlim([starttime,endtime]);
ylim([1,nTrial+1]);

%% labels
if isfiglabels
    title('','fontsize',fntsz);
    figset(gca,'time (sec)','trial index',fntsz);
else
    xlabel('');
    yyaxis left;
    ylabel('');
    yyaxis right;
    ylabel('');
    title('');
end

%% plot firing rate matrix
if isplotfr
    maxfr = max(lambda,[],'all');
    for i = 1:m
        subplot(3,m,2*m+i);
        b=barh(lambda(:,i));
        b.FaceColor = setcol(i,:);
        b.FaceAlpha = 0.5;
        b.EdgeColor = setcol(i,:);
        figset(gca,'','',fntsz);
        xlim([0,maxfr+5]);
        % label firing rate next to each bar
        for i1=1:N
            text(lambda(i1,i),i1,num2str(lambda(i1,i),'%0.1f'),...
                'HorizontalAlignment','left',...
                'VerticalAlignment','middle','FontSize',fntsz-2,'color',setcol(i,:))
        end
        if isfiglabels
            if i == ceil(m/2)
                xlabel('firing rate (Hz)');
            end
            if i == 1
                ylabel('neuron idx');
            end
        end
    end
end
end