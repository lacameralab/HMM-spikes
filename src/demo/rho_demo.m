% To quantify the variance captured by the models, we define D(HMM model) = sum_t=1^T
% sum_n=1^N (f_data(n,t)-dt*lambda_HMM(n,t))^2, summing over all time bins and
% neurons. f_data(n,t) is the spike count of neuron n in bin t,
% dt is the bin size, and lambda_HMM(n,t) is the estimated firing rate of neuron
% n in bin t. Then, we define \rho(HMM1, HMM2)=D(HMM1)/D(HMM2). If \rho
% is larger than 1, model HMM1 captures more variance than model HMM2, vice
% versa.
%
% This script shows how to calculated \rho for a trained HMM model (HMM1) and the
% true model (HMM2=HMM_true) using an MMPP (Markov-modulated Poisson
% process) dataset with 20 neurons and 10 states (50 trials, split into 40
% and 10 in the training and validation sets, respectively).
% 
% It also shows how to randomly shuffle the statese firing rates and the
% off-diagonal transition probabilities using the trained HMMs as examples.
% It will print out the the rho value of the trained HMM and that of their
% premuted model as a comparison.
% 
% This script also generates Figure 9C in the paper (the order of the
% panels may be different).
% 
% FUNCTIONALITY:
% 1. Loads MMPP dataset with known ground truth parameters
% 2. Calculates rho between trained models and true model
% 3. Performs state matching using Hungarian algorithm for visualization
% 4. Visualizes decoding and firing rate vectors for trained HMM and the true model
% 5. Compares rho value of the trained HMM against randomly permuted models
%
% Author: Tianshu Li
% Date: May 12, 2025

%% Environment Setup
% Add required paths
if ~exist('rasterPlot','file'), addpath '../utils'; end
if ~exist('phmm_stateProbs','file'), addpath '../pHMM'; end

clear;
close all;
% rng(12345); % Uncomment for reproducibility

%% Parameters
fntsz = 18; % font size for plots

% Set color scheme for state visualization
setcol = [
    0.1172, 0.5625, 0.9961;
    0.9961, 0.6445, 0.0000;  % orange
    0.0000, 0.0000, 1.0000;  % pure blue (kept for contrast)
    0.8000, 0.0000, 0.0000;  % red
    0.9020, 0.3800, 0.5880;  % pink
    0.7520, 0.8750, 0.0000;  % lime green
    0.1198, 0.3620, 0.2266;  % dark green
    0.6510, 0.3370, 0.1570;  % brown
    0.7194, 0.5208, 0.7194;  % lavender
    0.0000, 0.0000, 0.0000;  % black (anchor / neutral)
    0.6,0.6,0.6;
    1.0000, 0.2706, 0.0000;  % vivid orange-red (warm, distinct from red and orange)
    0.1216, 0.4706, 0.7059;  % teal-cyan (cool, distinct from blue and green)
];

col = lines;
setcol = [setcol; col];  % Add MATLAB's standard colors if needed

%% 1. Load the HMM models HMM1 and HMM2 (true model).
% Let's start with a MMPP dataset where we know the ground truth so that we
% can use the ground truth as the true model for comparison.

% 1.1 Load the MMPP dataset used to train HMM and the true model
HMMdataname = 'MMPP_N20m10T15_50trials_MCdt50ms_maxfr30_initstate0_3';
HMMdatadir = ['../../data/',HMMdataname,'.mat'];
HMMdata = load(HMMdatadir); % the MMPP dataset
% The true model is:
truemodel = struct();
truemodel.Lambda = HMMdata.MCpar.Lambda; % an N*m matrix, firing rate matrix, unit: Hz
truemodel.Gamma = HMMdata.MCpar.Gamma; % an m*m matrix, transition probability matrix

% Extract basic parameters
T = HMMdata.basicpar.T; % duration of each trial (P.S. this value is fixed for MMPP. However, it can be various for experimental datasets.)
nTrial = HMMdata.basicpar.nTrial; % total number of trials
dt = HMMdata.basicpar.dt; % the size of the time bin used for training HMM models
m = HMMdata.basicpar.m; % true number of states
N = HMMdata.basicpar.N; % number of neurons

% 1.2 Load some trained HMM model with the number of states 5, 10, and 12.
% Get the files containing the trained models
HMMmodeldir = '../../data/trained_HMM/';
f_HMMmodels = dir([HMMmodeldir,'*.mat']);
nHMMmodels = length(f_HMMmodels);

% Create figures
figure(1000); clf;

%% Visualize decoding by the true model
spikeshape = 'dot'; % the shape of spike, 'vertline' for vertical short line (raster), or 'dot'. 'dot' is faster

plottrial = 1; % the trial we want to plot

% Visualize decoding using the true model for selected trials
for k = plottrial
    % Calculate state probabilities
    pstates = phmm_stateProbs(HMMdata.spkc{k},m,truemodel.Lambda*dt,truemodel.Gamma);

    % Create figure for visualization
    figure(k); clf;
    subplot(nHMMmodels+1,1,1); cla;

    % Plot decoded states and spikes
    plotpHMM(HMMdata.MMPP(k).firings,pstates,N,HMMdata.pHMMtimev,'isnewfig',0,'fontsize',fntsz,'spikeshape',spikeshape,'statecolor',setcol);
    hold on;

    % Set figure properties
    figset(gca,'','',18);
    ylabel('neuron index','Interpreter','latex');
    yyaxis right;
    ylabel('state probability','Interpreter','latex');
    title(sprintf('true model, m = %i',m),'fontsize',fntsz,'Interpreter','latex');
end

%% Plot states (firing rate vectors, true model)
statealpha = 0.5;
stateidx = 1:m;
figure(1000);

% Plot each state's firing rate vectors
for jj = 1:length(stateidx)
    subplot(nHMMmodels+1,m,jj); cla;
    b=barh(truemodel.Lambda(:,stateidx(jj)));
    b.FaceColor = [0,0,0];
    b.EdgeColor = [0,0,0];
    figset(gca,'','',fntsz);
    set(gca,'color',[setcol(stateidx(jj),:),statealpha])
    ylim([0.5,N+0.5]);

    % Label only leftmost subplot
    if jj ~= 1
        set(gca,'YTickLabel',[]);
    else
        ylabel('neuron index','FontSize',fntsz);
    end
end

%% Process Each Trained HMM Model
for km = 1:nHMMmodels
    % Load the trained HMM model
    trainedHMM = load(fullfile(HMMmodeldir, f_HMMmodels(km).name));
    trainedmodel.Lambda = trainedHMM.estlambda; % an N*m matrix, firing rate matrix, unit: Hz
    trainedmodel.Gamma = trainedHMM.estGamma; % an m*m matrix, transition probability matrix
    M = trainedHMM.m; % the trained HMM's number of states

    %% Calculate \rho
    %% 1. Calculate the estimated firing rate in each time bin
    % To calculate D, we first calculate f_data(n,t), the spike count in each
    % time bin for each neuron.
    f_data = {}; % estimated firing rate by the data. Each cell is a trial.
    for ktrial = 1:nTrial
        f_data{ktrial} = HMMdata.spkc{ktrial}; % each row is one neuron, and each column is one time bin
    end

    %% 2. Find the state sequences decoded by the models.
    % To calulate D, we need to know which state each time bin is in and get
    % the firing rate of that state. For the model, we decode the data to get
    % the state sequences. For the ground truth, we take the ground truth
    % model as a trained HMM model and decode the state sequence to get
    % D(truemodel). We use the true model to decode the data instead of the
    % true state sequence to account for the variability due to state
    % transitions occurring inside a bin.

    isviterbi = false; % if true, use viterbi to decode; otherwise, assign the bin to the state with the largest posterior probability. Both methods provide the same result because their decoding is very similar.
    if isviterbi
        % Decode the data using trainedmodel to get the decoded state sequence by this model.
        % Viterbi decoding 
        trainedmodel_states = {}; % each cell is the state sequence of a trial. Use cell vector instead of a matrix because the length of each trial can vary.
        for ktrial = 1:nTrial
            trainedmodel_states{ktrial} = pHMM_viterbi(HMMdata.spkc{ktrial},M,trainedmodel.Lambda*dt,trainedmodel.Gamma);
        end

        % Decode the data using the truemodel to get the decoded state sequence by this model.
        % Viterbi decoding
        truemodel_states = {}; % each cell is the state sequence of a trial. Use cell vector instead of a matrix because the length of each trial can vary.
        for ktrial = 1:nTrial
            truemodel_states{ktrial} = pHMM_viterbi(HMMdata.spkc{ktrial},m,truemodel.Lambda*dt,truemodel.Gamma);
        end
    else % posterior decoding
        % Decode the data using trainedmodel to get the decoded state sequence by this model.
        trainedmodel_states = {}; % each cell is the state sequence of a trial. Use cell vector instead of a matrix because the length of each trial can vary.
        minp = 0; % only when the probability of the largest state is no less than minprob, the bin will be assigned to that state. Use 0 here to force a state to be assigned to all bins.
        for ktrial = 1:nTrial
            pstates = phmm_stateProbs(HMMdata.spkc{ktrial},M,trainedmodel.Lambda*dt,trainedmodel.Gamma);
            [pmax, maxstate] = max(pstates);
            trainedmodel_states{ktrial} = zeros(size(maxstate));
            trainedmodel_states{ktrial}(pmax>=minp) = maxstate(pmax>=minp);
        end
        % Decode the data using the truemodel to get the decoded state sequence by this model.
        truemodel_states = {}; % each cell is the state sequence of a trial. Use cell vector instead of a matrix because the length of each trial can vary.
        for ktrial = 1:nTrial
            pstates = phmm_stateProbs(HMMdata.spkc{ktrial},m,truemodel.Lambda*dt,truemodel.Gamma);
            [pmax, maxstate] = max(pstates);
            truemodel_states{ktrial} = zeros(size(maxstate));
            truemodel_states{ktrial}(pmax>=minp) = maxstate(pmax>=minp);
        end
    end

    %% 3. Calculate the variance captured by the model on the validation set (D)
    % D(HMM model) = sum_t=1^T sum_n=1^N (f_data(n,t)-dt*lambda_HMM(n,t))^2
    trialidx = trainedHMM.val_ind; % Validation set. Use trainedHMM.train_ind for training set

    % 4.1 the trained HMM model (trainedmodel)
    D_trainedmodel = 0;
    for ktrial = trialidx
        D_trainedmodel = D_trainedmodel + sum((f_data{ktrial} - dt*trainedmodel.Lambda(:,trainedmodel_states{ktrial})).^2,'all');
    end

    % 4.2 the true model (truemodel)
    D_truemodel = 0;
    for ktrial = trialidx
        D_truemodel = D_truemodel + sum((f_data{ktrial} - dt*truemodel.Lambda(:,truemodel_states{ktrial})).^2,'all'); 
    end

    % \rho
    rho_HMMvsTrue = D_trainedmodel/D_truemodel;


    %% Match states using Hungarian algorithm for Visualization
    newIndices = optimal_pairing(trainedmodel.Lambda, truemodel.Lambda); % paring of the states in two models, the first column is HMMmodel1, the second column is truemodel.
    [a, I] = sort(newIndices(:,1)); % the new indices of the states in HMMmodel1 w.r.t. the ones in truemodel.
    
    % Map trained model states to true model states for consistent colors
    coloridx = nan+zeros(1,M);
    coloridx(a) = newIndices(I,2)';
    if length(a) < M
        coloridx(isnan(coloridx)) = m+1:M;
    end

    %% =======================     Visualization    ===========================
    %% Decoding (trained HMM)
    for k = plottrial
        pstates = phmm_stateProbs(HMMdata.spkc{k},M,trainedmodel.Lambda*dt,trainedmodel.Gamma);

        figure(k);
        subplot(nHMMmodels+1,1,km+1); cla;
        plotpHMM(HMMdata.MMPP(k).firings,pstates,N,HMMdata.pHMMtimev,'isnewfig',0,'fontsize',fntsz,'spikeshape',spikeshape,'statecolor',setcol(coloridx,:));
        hold on;

        % Set figure properties
        figset(gca,'','neuron indx',18);
        if km == nHMMmodels
            xlabel('time (sec)','Interpreter','latex');
        end
        ylabel('neuron index','Interpreter','latex');
        yyaxis right;
        ylabel('state probability','Interpreter','latex');
        title(sprintf('$m = %i, \\rho = %.3f$',M,rho_HMMvsTrue),'fontsize',fntsz,'Interpreter','latex');
    end

    %% States (trained HMM)
    statealpha = 0.5;
    stateidx = 1:M;
    figure(1000);

    % Plot each state's firing rate vectors
    for jj = 1:length(stateidx)
        subplot(nHMMmodels+1,M,km*M+jj); cla;
        b=barh(trainedmodel.Lambda(:,stateidx(jj)));
        b.FaceColor = [0,0,0];
        b.EdgeColor = [0,0,0];
        figset(gca,'','',fntsz);
        set(gca,'color',[setcol(coloridx(stateidx(jj)),:),statealpha])
        ylim([0.5,N+0.5]);

        % Add labels to appropriate subplots
        if jj ~= 1
            set(gca,'YTickLabel',[]);
        else
            ylabel('neuron index','FontSize',fntsz);
        end
        if km == nHMMmodels && jj == 1
            xlabel('firing rate (Hz)','FontSize',fntsz);
        end
    end

    %% =================================================================
    % For comparison, we randomly shuffling the states’ firing rates and the off-diagonal
    % transition probabilities. This is the value expected by chance when comparing the true
    % model with a model with the same firing rates and transition rates as
    % the test model.

    %% Randomly shuffle Lambda and Gamma
    permuted_HMM = permuteHMMmodel(trainedmodel);
    if isviterbi
        % Decode the data using permuted_HMM to get the decoded state sequence by this model.
        % Viterbi decoding
        permuted_HMM_states = {}; % each cell is the state sequence of a trial. Use cell vector instead of a matrix because the length of each trial can vary.
        for ktrial = 1:nTrial
            permuted_HMM_states{ktrial} = pHMM_viterbi(HMMdata.spkc{ktrial},size(permuted_HMM.Lambda,2),permuted_HMM.Lambda*dt,permuted_HMM.Gamma);
        end
    else % posterior decoding
        % Decode the data using HMMmodel1 to get the decoded state sequence by this model.
        permuted_HMM_states = {}; % each cell is the state sequence of a trial. Use cell vector instead of a matrix because the length of each trial can vary.
        minp = 0; % only when the probability of the largest state is no less than minprob, the bin will be assigned to that state. Use 0 here to force a state to be assigned to all bins.
        for ktrial = 1:nTrial
            pstates = phmm_stateProbs(HMMdata.spkc{ktrial},M,permuted_HMM.Lambda*dt,permuted_HMM.Gamma);
            [pmax, maxstate] = max(pstates);
            permuted_HMM_states{ktrial} = zeros(size(maxstate));
            permuted_HMM_states{ktrial}(pmax>=minp) = maxstate(pmax>=minp);
        end
    end

    % the permuted model (permuted_HMM)
    D_HMMmodel_perm = 0;
    for ktrial = trialidx
        D_HMMmodel_perm = D_HMMmodel_perm + sum((f_data{ktrial} - dt*permuted_HMM.Lambda(:,permuted_HMM_states{ktrial})).^2,'all');
    end
    rho_HMMvsTrue_perm = D_HMMmodel_perm/D_truemodel;

    % Print out the rho value of the trained HMM and the permuted HMM as a comparison
    fprintf('%i states: rho(trained HMM) = %.3f, rho(permuted HMM) = %.3f\n',M,rho_HMMvsTrue,rho_HMMvsTrue_perm);

end
