function [stateseq] = pHMM_viterbi(x,m,lambda ,Gamma ,delta)

if nargin <5
    delta = ones(1,m)/(eye(m)-Gamma+ones(m)); % stationary distribution
end
n = length(x);
P_xt = zeros(m,n); % P(xt), p(i,t) = P(xt|Ci)
for i = 1:m
    P_xt(i,:) = exp(sum(-lambda(:,i)+log(lambda(:,i)).*x-log(factorial(x)),1));
end
xi = zeros(n,m);
xi1 = delta*P_xt(:,1);
xi(1,:) = xi1/sum(xi1);

for t = 2:n
%     xit = max(xi(t-1,:)*Gamma)*P_xt(:,t);
%     xi(t,:) = xit/sum(xit);
    xi_t = zeros(1,m);
    for j = 1:m
        xi_t(j) = max(xi(t-1,:).*Gamma(:,j)')*P_xt(j,t);
    end
    xi(t,:) = xi_t/sum(xi_t);
end
stateseq = zeros(1,n);
[x,stateseq(n)] = max(xi(n,:));
for t = (n-1):-1:1
    [x, stateseq(t)] = max(Gamma(:,stateseq(t+1)).*xi(t,:)');
end
end