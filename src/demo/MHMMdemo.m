% This file is a demo for applying Multinoulli-HMM to spiking data.
% 
% NOTICE:
% The unit of time is SECOND though out the script. Please change the unit
% of time in your data to second before running this script!
% 
% Author: Tianshu Li
% Date: 10/25/2021

rng(3456); % for reproducibility, remove it to get different initial conditions

clear;
addpath '../utils';
addpath '../mHMM';

%  -------------   Set color   ---------------
cols = lines;
setcol = cols;
% set color for each state
setcol(1:8,:) = ...
         [0.9961    0.6445         0
          0.8500    0.3250    0.0980
          0         0         1.0000
          0.1198    0.3620    0.2266
          0.1172    0.5625    0.9961
          0.8000         0         0
          0.7520    0.8750         0
          0.7194    0.5208    0.7194];
% --------------------------------------------

ischecking = 0; % if 1, enables code checking warnings

if ischecking
    prompt = '======\nThe unit of time is SECOND though out the script. \nIs every time unit in your data AND the script (after your edit) SECOND? \n(type ''yes''/''no'')\n';
    istimeunitcorrect = input(prompt);
    if strcmpi(istimeunitcorrect,'no')
        error('Please change the unit of time in your data AND the script (after your edit) to SECOND before running this script!');
    end
end

%% =======================      Load data     =============================
% TO DO:
% Put the directory of your spike data here
% Requirement:
%   1. time unit should be second.
%   2. the firing data should be stored in a structure "data",
%      data(k).firings is a two column matrix storing all the recorded
%      neurons' spike time in trial k. The first column is the spike time
%      and the second column is the index of the neuron which fires at the
%      corresponding time on the same row.
%   3. variable "starttime" is a scalar which is the starting time of a
%      trial (the same for all trials). 
%   4. variable "endtime" is a scalar which is the end time of a
%      trial (the same for all trials). 
%   5. variable "N" is a scalar which is the number of neurons recorded in
%      data.

datadir = '../../data/exampledata.mat'; % * put the directory of your prepared spiking data here *
datafilename = split(datadir,'/');
datafilename = datafilename{end}; % name of data file
load(datadir);
% data=data(1:10); % select only first 10 trials

%% ===========    Seperate training set and validation set    =============
% TO DO:
% 1. You can change the percentage/number of trials for training set (and
%    validation set) by changing 'ptrain'/'ntrain'.
% 2. You can set the protocol of selecting trials for training set (and
%    validation set) by changing 'train_ind'.

ptrain = 0.8; % * the percentage of trials used as training set *
nTrial = length(data); % total number of trials
ntrain = ceil(nTrial*ptrain); % number of trials in training set
nval = nTrial-ntrain; % number of trials in validation set

% **********
% Choosing trials for training set. This demo randomly chooses 'ntrain'
% trials to form training set. You can change the way of choosing trials
% for training set, especially when the experimental setting of your
% dataset cannot be assumed as identical for all trials (e.g. different
% external stimuli were provided). 
% **********
train_ind = sort(randperm(nTrial,ntrain)); % trial indices in training set
val_ind = setdiff(1:nTrial,train_ind); % trial indices in validation set

%% =======    Generate observation sequences for Multinoulli-HMM    =======
% TO DO:
%   1. You can change the time bin 'dt' for each observation (a vector of
%   spike counts for each neuron in that time bin).
%   2. You can choose the neurons used for Poisson-HMM by changing 'neuind'.

dt = 0.005; % * sec, size of time bin, suggest 0.005 sec *
timev = starttime:dt:endtime; % time vector

neuind = 1:N; % selected neurons' indices
nNeuron = length(neuind); % number of neurons selected for Poisson-HMM

seq = []; % observation sequences for Multinoulli-HMM, each row is a trial, each entry in a row is the index of the neuron which fires at the corresponding time bin (if no neuron fires, the value is nNeuron+1)
for i = 1:nTrial
    seq = [seq; firings2seq(data(i).firings,neuind,timev)];
end

%% =========    Choose initial conditions for Multinoulli-HMM    ==========
% TO DO:
%   1. You can decide the number of states by changing 'm'.
%   2. You can change the way of choosing the initial parameters:
%      transition matrix (TRANS_GUESS) and emission probability matrix
%      (EMIS_GUESS). Multinoulli-HMM strongly depends on the initial
%      conditions. This demo provided a good pair.

m = 10; % number of states

TRANS_GUESS = genTRANS(m,1);  % initial transition matrix
EMIS_GUESS = 0.1*rand(m,N+1)/N; EMIS_GUESS(:,1)=1-sum(EMIS_GUESS(:,2:end),2); % initial emission probability matrix

%% ======================     Multinoulli-HMM     =========================
[estTRANS,estEMIS,converged,logliks,rec] = myhmmtrain(seq(train_ind,:),TRANS_GUESS,EMIS_GUESS,'maxiterations',1000,'tolerance',1e-4);

estTRANS
estEMIS

%% =======================     Visualization    ===========================
fntsz = 20; % fontsize

%% self transition prob
figure(1000); clf; hold all;
for i = 1:m
    plot(rec.TR_diag(i,rec.TR_diag(i,:)>0));
end
legend(strcat('State ',string(1:m)),'location','best');
legend 'boxoff';
figset(gca,'iteration','self transition probability',fntsz);
title('Self transition probability during training (Multinoulli-HMM)');
hold off;

%% log-likelihood during training
figure(1001); clf;
plot(logliks);
figset(gca,'number of iterations','log-likelihood');
title('Multinoulli-HMM');

%% stationary distribution
disp('Stationary distribution:');
delta = ones(1,m)/(eye(m)-estTRANS+ones(m)) % stationary distribution

%% decoding
spikeshape = 'vertline'; % the shape of spike, 'vertline' for vertical short line (raster), or 'dot'. 'dot' is faster

for k = 1:min(3,nTrial)
    figure(k); clf; hold on;
    set(gcf,'position',[1681          99        1292         725]);
    subplot(2,1,1); cla;
    plotmHMM(seq(k,:), estTRANS, estEMIS, data(k).firings,timev,'isnewfig',0,'statecolor',setcol,'spikeshape',spikeshape);
    figset(gca,'time (sec)','neuron indx',fntsz);
    title(sprintf('Multinoulli-HMM (N=%i), m=%i, dt=%i ms, trial %i',N,m,dt*1000,k),'fontsize',fntsz+2,'fontweight','bold');
    hold off;
    
    % viterbi decoding
    stateseq = hmmviterbi(seq(k,:),estTRANS,estEMIS);
    subplot(2,1,2); cla;
    plotviterbi(data(k).firings,stateseq,nNeuron,timev,'isnewfig',0,'statecolor',setcol,'spikeshape',spikeshape);
    figset(gca,'time (sec)','neuron indx',18);
    title(sprintf('Viterbi decoding'));
end
%% plot all trials in one figure
isseptrainvalset = 1; % if 1, plot all training trials before validation trials; if 0, plot selected trials as the order in 'sidx'
sidx = 1:nTrial; % choose the trials to be plotted and the order of them
stateshape = 'square';%'probcurve'; % the presentation of states. If 'square', each state is colored as a square between the trial indices. If 'probcurve', only the area under the probability curve of the estimated state will be colored.

figure(1002); clf; hold on;
if isseptrainvalset
    plotmHMMtrials(seq(sort([train_ind, val_ind]),:),estTRANS,estEMIS,timev,'isnewfig',0,'fontsize',fntsz,'statecolor',setcol,'stateshape',stateshape);
    plot([starttime,endtime], 1+[ntrain,ntrain],'k-','linewidth',2); % seperate training set and validation set
else
    plotmHMMtrials(seq(sidx,:),estTRANS,estEMIS,timev,'isnewfig',0,'fontsize',fntsz,'statecolor',setcol,'stateshape',stateshape);
end
title(sprintf('Multinoulli-HMM, %s',replace(datafilename,'_','\_')));