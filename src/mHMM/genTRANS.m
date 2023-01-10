function [TRANS] = genTRANS(nState,isbiased,a)
% [TRANS] = genTRANS(nState,isbiased,a)
%   This function generates a random transition probability matrix.
% Inputs:
%   nState - number of state.
%   isbiased - bool, if isbiased == 1, generate TRANS whose diagonal values
%              are close to one. Default: 0.
%   a - double, a scaling parameter which controls the diagonal of
%       transition probability matrix. (The larger a is, the lower the
%       diagonal will be.) Default: 0.01.
% Outputs:
%   TRANS - symmetrical nState*nState matrix, TRANS(i,j)=P(state i->state
%   j), TRANS(i,j)==TRANS(j,i), TRANS(i,i)~=0, the sum of each row is 1.
% 
% - Tianshu Li, 8/9/2019

if nargin < 2
    isbiased = 0;
end
if nargin < 3
    a = 0.01;
end
if isbiased
    TRANS = eye(nState)+a/nState*rand(nState);
else
    TRANS = rand(nState);
end
m = size(TRANS,1);

for i = 1:nState
    % normalize the first row(i) that has not been normalized, then
    % symmitrize the small matrix start at (i,i).
    TRANS(i,i:end) = (1-sum(TRANS(i,1:(i-1))))*TRANS(i,i:end)/sum(TRANS(i,i:end));
    TRANS = triu(TRANS)+triu(TRANS)'-eye(m).*diag(TRANS);
end

end
