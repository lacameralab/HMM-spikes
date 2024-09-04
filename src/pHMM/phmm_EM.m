function [lambda,Gamma,isconverge,deltas,logliks,rec,iter,crit,llk,AIC,BIC] = phmm_EM(xs,m,lambda,Gamma,thresh,maxiter,tol)
% [lambda,Gamma,isconverge,deltas,logliks,rec,iter,crit,llk,AIC,BIC] = phmm_EM(xs,m,lambda,Gamma,thresh,maxiter,tol)
% Sticky Poisson HMM of spike data (set 'thresh' to zero for a 'standard' Poisson-HMM).
%
% Inputs:
%   xs: observation vector, if there are several trials, put them in a
%       cell.
%   m: number of states
%   lambda: firing rate (parameter of Poisson distribution) in each state
%   Gamma: transition probability matrix (Gamma(i,j)=P(State_i-->State_j))
%   thresh: threshold of the self-transition probability (diagonal of
%       Gamma). Set thresh to zero for a 'standard' Poisson-HMM
%   maxiter: maxiteration number. Default: 1000
%   tol: tolerence. Default: 1e-6
% Outputs:
%   lambda: estimated firing rate matrix (lambda/dt will have unit Hz).
%   Gamma: estimated transition probability matrix.
%   isconverge: 1 if the algorithm is converged, otherwise 0.
%   deltas: each row is the estimated initial distribution of the
%       correspongding trial.
%   loglikes: a vector contains log-likelihood of each iteration.
%   rec: recordings, diagonal of Gamma across iteration and firing rate of
%       two neurons across iteration.
%   iter: iteration reached when return.
%   crit: convergence criterion, Euclid distance (norm 2)
%       distance(oldGamma,Gamma)+distance(oldlambda,lambda).
%   llk: a scalar, log-likelihood of the estimated model on the training set.
%   AIC: a scalar, Akaike information criterion
%   BIC: a scalar, Bayesian information criterion
% 
% Tianshu Li - inspired by the R code in Zucchini's book A.2.3.

if nargin < 7
    tol = 1e-4; % tolerence
end
if nargin < 6
    maxiter = 1000; % max number of iteration
end

if nargin < 5
    thresh = 0.8;
end

if ~isempty(find(diag(Gamma)<thresh,1))
    error('The diagonal of transition matrix should be larger than thresh (default: 0.8)!');
end

if ~iscell(xs)
    xs = {xs};
end

ntrial = length(xs);
deltas = ones(ntrial,m)/m;
N = size(lambda,1);
Gamma_init = Gamma;

rec.Gamma_diag = zeros(m,maxiter*ntrial);
rec.lambda_neu1 = zeros(m,maxiter*ntrial);
l = 1; % track the step of update
rec.Gamma_diag(:,l) = diag(Gamma);
rec.lambda_neu1(:,l) = lambda(1,:)';

flag = 0;
tauref = 50; % iterations of the refractory period
isconverge = 0;
ref = tauref; % refractory period after each hit reset
rec.hit = []; % record the iteration and trial where the threshold was hit
latestgoodGamma = Gamma; % the latest Gamma whose diagonal is above the threshold
latestgoodlambda = lambda; % the corresponding lamdba of latestgoodGamma

loglik = 1;
logliks = zeros(1,maxiter);

% To prevent repeatly calculating log(factorial(x)), store the result ahead of time.
LFxs = cell(1,ntrial);
for tl = 1:ntrial
    LFxs{tl} = log(factorial(xs{tl}));
end

for iter = 1:maxiter
    oldLL = loglik;
    oldlambda = lambda;
    oldGamma = Gamma;
    loglik = 0;
    Gamma_next = zeros(m);
    lambda_next = zeros(N,m);
    lambdascale = zeros(1,m);
    Gammascale = zeros(m,1);
    ref = ref - 1;
    l = l+1;
    for tl = 1:ntrial
        x = xs{tl};
        n = size(x,2);	
        delta = deltas(tl,:);
        
        % P(xt)
        lP_xt = zeros(m,n); % P(xt), p(i,t) = P(xt|Ci)
        for i = 1:m
            lP_xt(i,:) = sum(-lambda(:,i)+log(lambda(:,i)).*x-LFxs{tl},1);
        end

        [la, lb] = phmm_lalphabeta(m,x,oldlambda,oldGamma,delta,lP_xt);
        c = max(la(:,end));
        llk = c + log(sum(exp(la(:,end)-c)));
        loglik = loglik + llk;
        for j = 1:m
            for k = 1:m
                Gamma_next(j,k) = Gamma_next(j,k) + oldGamma(j,k)*sum(exp(la(j,1:n-1)+lP_xt(k,2:n)+lb(k,2:n)-llk));
            end
            lambda_next(:,j) = lambda_next(:,j) + sum(exp(la(j,:)+lb(j,:)-llk).*x,2);
            lambdascale(j) = lambdascale(j) + sum(exp(la(j,:)+lb(j,:)-llk));
        end
        
        Gammascale = sum(Gamma_next,2);
        delta_next = exp(la(:,1)+lb(:,1)-llk);
        delta_next = delta_next'/sum(delta_next);
        
        deltas(tl,:) = delta_next;
        
    end
    
    Gamma_next = Gamma_next./Gammascale;
    Gamma_next = Gamma_next./sum(Gamma_next,2);
    lambda_next = lambda_next./lambdascale;
    
    % impose the constrain
    indG = find(diag(Gamma_next)<thresh);
    if 1 && ((~isempty(indG))) && all(abs(diag(Gamma_next)-diag(oldGamma))<5e-5,'all')
        Gamma_next = latestgoodGamma; % Gamma
        lambda_next = latestgoodlambda; % lambda
        lambda_next(:,indG) = latestgoodlambda(randperm(N),indG);
        flag = 1;
        ref = tauref; % refractory period after each hit reset
        rec.hit = [rec.hit; [iter,tl]];
    end
    indL = find(lambda_next==0);
    if 1 && ~isempty(indL)
        lambda_next(indL) = 0.001;
    end
    if 1 && isempty(indG)
        latestgoodGamma = Gamma_next;
        latestgoodlambda = lambda_next;
    end
        
    % record the update of diagonal of transition matrix
    rec.Gamma_diag(:,l) = diag(Gamma_next);
    % record the update of the firing rate of neuron 1 and 2
    rec.lambda_neu1(:,l) = lambda_next(1,:)';
    rec.lambda_neu2(:,l) = lambda_next(2,:)';
    
    % crit: convergence criterion (can be replaced by some other ceiterion
    % chosen by user).
    crit = norm(oldlambda-lambda_next) + norm(oldGamma-Gamma_next);
    critd = 0;
    crit = crit + critd + flag;
    flag = 0;
    
    logliks(iter) = loglik;
    if abs(loglik-oldLL) < 1e-6
        if crit < tol
            if isempty(indG)
                np = m*m+m-1;
                AIC = -2*(llk-np);
                BIC = -2*llk+np*log(n);
                isconverge = 1;
                return;
            else
                ref = 0;
            end
        end
    end
    
    lambda = lambda_next; %lambda_next;
    Gamma = Gamma_next; %Gamma_next;
end

np = m*m+m-1;
AIC = -2*(llk-np);
BIC = -2*llk+np*log(n);
fprintf('No convergence after %i iteration.\n', maxiter);

if Gamma == Gamma_init
    fprintf('No change of Gamma after %i iteration.\n', maxiter);
end

end
