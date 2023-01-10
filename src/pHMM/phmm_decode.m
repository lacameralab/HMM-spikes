function [pStates,pSeq, fs, bs, s] = phmm_decode(x,Gamma,lambda,varargin)
% phmm_decode is the Poisson HMM version of MATLAB function hmmdecode.
% Gamma is tr in hmmdecode, P_xt is e in hmmdecode.
% Gamma must be square
% 
% HMMDECODE calculates the posterior state probabilities of a sequence.
%   PSTATES = HMMDECODE(SEQ,TRANSITIONS,EMISSIONS) calculates the posterior
%   state probabilities, PSTATES, of sequence SEQ from a Hidden Markov
%   Model specified by transition probability matrix,  TRANSITIONS, and
%   EMISSION probability matrix, EMISSIONS. TRANSITIONS(I,J) is the
%   probability of transition from state I to state J. EMISSIONS(K,SYM) is
%   the probability that symbol SYM is emitted from state K. The posterior
%   probability of sequence SEQ is the probability P(state at step i = k |
%   SEQ).  PSTATES is an array with the same length as SEQ and one row for
%   each state in the model. The (i,j) element of PSTATES gives the
%   probability that the model was in state i at the jth step of SEQ.
%
%   [PSTATES, LOGPSEQ] = HMMDECODE(SEQ,TR,E) returns, LOGPSEQ, the log of
%   the probability of sequence SEQ given transition matrix, TR and
%   emission matrix, E.
%
%   [PSTATES, LOGPSEQ, FORWARD, BACKWARD, S] = HMMDECODE(SEQ,TR,E) returns
%   the forward and backward probabilities of the sequence scaled by S. The
%   actual forward probabilities can be recovered by using:
%        f = FORWARD.*repmat(cumprod(s),size(FORWARD,1),1);
%   The actual backward probabilities can be recovered by using:
%       bscale = cumprod(S, 'reverse'); b =
%       BACKWARD.*repmat([bscale(2:end), 1],size(BACKWARD,1),1);
%
%   HMMDECODE(...,'SYMBOLS',SYMBOLS) allows you to specify the symbols that
%   are emitted. SYMBOLS can be a numeric array or a string array cell
%   array of the names of the symbols.  The default symbols are integers 1
%   through M, where N is the number of possible emissions.
%
%   This function always starts the model in state 1 and then makes a
%   transition to the first step using the probabilities in the first row
%   of the transition matrix. So in the example given below, the first
%   element of the output states will be 1 with probability 0.95 and 2 with
%   probability .05.
%
%   Examples:
%
%       tr = [0.95,0.05; ...
%             0.10,0.90];
%           
%       e = [1/6,  1/6,  1/6,  1/6,  1/6,  1/6; ...
%            1/10, 1/10, 1/10, 1/10, 1/10, 1/2;];
%
%       [seq, states] = hmmgenerate(100,tr,e); 
%       pStates = hmmdecode(seq,tr,e);
%
%       [seq, states] = hmmgenerate(100,tr,e,'Symbols',...
%                 {'one','two','three','four','five','six'});
%       pStates = hmmdecode(seq,tr,e,'Symbols',...
%                 {'one','two','three','four','five','six'});
%
%   See also HMMGENERATE, HMMESTIMATE, HMMVITERBI, HMMTRAIN.
% 
%   Reference: Biological Sequence Analysis, Durbin, Eddy, Krogh, and
%   Mitchison, Cambridge University Press, 1998.
% 
%   Copyright 1993-2014 The MathWorks, Inc.

if nargin > 0
    if isstring(x)
        x = cellstr(x);
    end
end

if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

numStates = size(Gamma,1);
checkTr = size(Gamma,2);
if checkTr ~= numStates
    error(message('stats:hmmdecode:BadTransitions'));
end

% number of columns of lambda must be same as number of states

checklambda  = size(lambda,2);
if checklambda ~= numStates
    error(message('stats:hmmdecode:InputSizeMismatch'));
end

% add extra data to start to make algorithm cleaner at f0 and b0
x = [zeros(size(lambda,1),1), x ];
T = size(x,2);

% P(xt) whis is e in hmmdecode
P_xt = zeros(numStates,T); % P(xt), p(i,t) = P(xt|Ci)
for i = 1:numStates
    P_xt(i,:) = exp(sum(-lambda(:,i)+log(lambda(:,i)).*x-log(factorial(x)),1));
%     P_xt(i,:) = exp(sum(-lambda(:,i)+log(lambda(:,i)).*x)); % it is not necessary to include the factorial term in the observation density because it is common to all states.
end

% This is what we'd like to do but it is numerically unstable
% warnState = warning('off');
% logTR = log(tr);
% logE = log(e);
% warning(warnState);
% f = zeros(numStates,L);
% f(1,1) = 1;
% % for count = 2:L
%     for state = 1:numStates
%         f(state,count) = logE(state,seq(count)) + log(sum( exp(f(:,count-1) + logTR(:,state))));
%     end
% end
% f = exp(f);

% so we introduce a scaling factor
fs = zeros(numStates,T);
fs(1,1) = 1;  % assume that we start in state 1.
s = zeros(1,T);
s(1) = 1;
for count = 2:T
    for state = 1:numStates
        fs(state,count) = P_xt(state,count) .* (sum(fs(:,count-1) .*Gamma(:,state)));
    end
    % scale factor normalizes sum(fs,count) to be 1. 
    s(count) =  sum(fs(:,count));
    fs(:,count) =  fs(:,count)./s(count);
end

%  The  actual forward and  probabilities can be recovered by using
%   f = fs.*repmat(cumprod(s),size(fs,1),1);


% This is what we'd like to do but it is numerically unstable
% b = zeros(numStates,L);
% for count = L-1:-1:1
%     for state = 1:numStates
%         b(state,count) = log(sum(exp(logTR(state,:)' + logE(:,seq(count+1)) + b(:,count+1)  )));
%     end
% end

% so once again use the scale factor
bs = ones(numStates,T);
for count = T-1:-1:1
    for state = 1:numStates
      bs(state,count) = (1/s(count+1)) * sum( Gamma(state,:)'.* bs(:,count+1) .* P_xt(:,count+1)); 
    end
end

%  The  actual backward and  probabilities can be recovered by using
%  scales = cumprod(s, 'reverse'); 
%  b = bs.*repmat([scales(2:end), 1],size(bs,1),1);

pSeq = sum(log(s));
pStates = fs.*bs;

% get rid of the column that we stuck in to deal with the f0 and b0 
pStates(:,1) = [];


