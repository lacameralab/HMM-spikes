% This file is a demo for applying sticky Poisson-HMM to spiking data.
% 
% NOTICE:
% The unit of time is SECOND throughout the script. Please change the unit
% of time in your data to second before running this script!
% 
% Author: Tianshu Li
% Date: 9/21/2021

rng(3456); % for reproducibility, remove it to get different initial conditions

clear;
addpath '../basic_functions';
addpath '../pHMM';

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
        error('Please change the units of time in your data AND the script (after your edit) to SECONDS before running this script!');
    end
end

%% =======================      Load data     =============================
% Put the directory of your spike data here
% Requirements:
%   1. time units should be seconds
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

%% =========    Generate observation sequences for Poisson-HMM    =========
%   1. You can change the time bin 'dt' for each observation (a vector of
%   spike counts for each neuron in that time bin).
%   2. You can choose the neurons used for Poisson-HMM by changing 'neuind'.

dt = 0.05; % * sec, size of time bin, suggest 0.05 sec *
timev = starttime:dt:endtime; % time vector

neuind = 1:N; % selected neurons' indices
nNeuron = length(neuind); % number of neurons selected for Poisson-HMM

spkc = {}; % observation sequences for Poisson-HMM, each entry is a matrix of spike count of every neuron is every time bin in one trial
for i = 1:nTrial
    if isempty(data(i).firings)
        spkc{i} = zeros(N,length(timev)-1);
    else
        spkc{i} = spikecount(data(i).firings, dt, timev, neuind);
    end
end


%% =============    Choose initial conditions for Poisson-HMM    ==========
%   1. You can decide the number of states by changing 'm'.
%   2. You can change the threshold of self-transition probability (detail
%      see paper) by changing 'thresh' (0.8 was used in paper).
%   3. You can change the way of choosing the initial firing rate matrix (lambda_guess).
%      This demo set the initial firing rate of each neuron in each state
%      to be random values between its maximum and minimum firing rate. For
%      example, you can also determine initial firing rate matrix by
%      clustering the firing rate vectors at all time bins.
%   4. You can change the way of setting the initial transition matrix (Gamma_guess).
%      This demo set the initial self-transition probability of each state
%      to be higher than 'thresh'. Each entry of Gamma_guess should be value
%      between 0 and 1. The sum of each row should be 1. 

m = 10; % number of states
thresh = 0.8; % threshold of self-transition probability (detail see paper)

% Find the maximum firing rate of each neuron in all trials. Use this value
% to generate initial firing rate matrix 
maxfr = zeros(1,nNeuron);
for i = 1:nNeuron
    for j = 1:nTrial
        maxfr(i) = max(maxfr(neuind(i)), max(firingrate(data(j).firings, neuind(i))));
    end
end

lambda_guess = maxfr'.*rand(nNeuron,m)*dt; % initial value for firing rate matrix
Gamma_guess = rand(m); Gamma_guess = thresh*eye(m)+(1-thresh)*Gamma_guess./sum(Gamma_guess,2); % initial value for transition matrix

%% ========================     Poisson-HMM     ===========================
[estlambda,estGamma,isconverge,delta1,logliks,rec,iter,crit] = phmm_EM(spkc(train_ind),m,lambda_guess,Gamma_guess,thresh);

estGamma
estlambda = estlambda/dt

%% =======================     Visualization    ===========================
fntsz = 20; % fontsize

%% self transition prob
figure(1000); clf; hold all;
for i = 1:m
    plot(rec.Gamma_diag(i,rec.Gamma_diag(i,:)>0),'color',setcol(i,:));
end
legend(strcat('State ',string(1:m)),'location','best');
legend 'boxoff';
figset(gca,'iteration','self transition probability',fntsz);
title('Self transition probability during training (sticky Poisson-HMM)');
hold off;

%% log-likelihood during training
figure(1001); clf;
plot(logliks(logliks~=0));
figset(gca,'number of iterations','log-likelihood');
title('sticky Poisson HMM');

%% decoding
spikeshape = 'vertline'; %'dot';% the shape of spike, 'vertline' for vertical short line (raster), or 'dot'. 'dot' is faster

for k = 1:min(3,nTrial)
    pstates = phmm_stateProbs(spkc{k},m,estlambda*dt,estGamma);
    
    figure(k); clf; hold on; 
    set(gcf,'position',[1681          99        1292         725]);
    subplot(2,1,1); cla;
    plotpHMM(data(k).firings,pstates,nNeuron,timev,'isnewfig',0,'fontsize',fntsz,'spikeshape',spikeshape,'statecolor',setcol);
    hold on;
    figset(gca,'time (sec)','neuron indx',18);
    title(sprintf('sticky Poisson-HMM (N=%i), m=%i, dt=%i ms, trial %i',N,m,dt*1000,k),'fontsize',fntsz+2,'fontweight','bold');
    
    % viterbi decoding
    stateseq = pHMM_viterbi(spkc{k},m,estlambda*dt,estGamma);
    subplot(2,1,2); cla;
    plotviterbi(data(k).firings,stateseq,nNeuron,timev,'isnewfig',0,'statecolor',setcol,'spikeshape',spikeshape);
    figset(gca,'time (sec)','neuron indx',18);
    title(sprintf('Viterbi decoding'));
end

%% stationary distribution
disp('Stationary distribution:');
delta = ones(1,m)/(eye(m)-estGamma+ones(m)) % stationary distribution

%% plot all trials in one figure
isseptrainvalset = 1; % if 1, plot all training trials before validation trials; if 0, plot selected trials as the order in 'sidx'
sidx = 1:nTrial; % choose the trials to be plotted and the order of them
stateshape = 'square';%'probcurve'; % the presentation of states. If 'square', each state is colored as a square between the trial indices. If 'probcurve', only the area under the probability curve of the estimated state will be colored.

figure(1002); clf; hold on;
if isseptrainvalset
    ah = plotpHMMtrials(spkc(sort([train_ind, val_ind])),estlambda,estGamma,timev,'isnewfig',0,'fontsize',fntsz,'statecolor',setcol,'stateshape',stateshape);
    plot(ah, [starttime,endtime], 1+[ntrain,ntrain],'k-','linewidth',2); % seperate training set and validation set
else
    plotpHMMtrials(spkc(sidx),estlambda,estGamma,timev,'isnewfig',0,'fontsize',fntsz,'statecolor',setcol,'stateshape',stateshape);
end
sgtitle(sprintf('sticky Poisson-HMM, threshold=%g, %s',thresh,replace(datafilename,'_','\_')));
