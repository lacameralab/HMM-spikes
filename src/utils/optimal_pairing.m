function assignment = optimal_pairing(A, B)
%OPTIMAL_PAIRING Finds the optimal pairing between two sets of vectors to minimize Euclidean distance
%
% USAGE:
%   assignment = optimal_pairing(A, B)
%
% INPUTS:
%   A - NxM1 matrix where each column represents an N-dimensional vector
%   B - NxM2 matrix where each column represents an N-dimensional vector
%
% OUTPUTS:
%   assignment - Mx2 matrix where each row [i,j] indicates that vector A(:,i)
%                is paired with vector B(:,j) in the optimal assignment
%
% DESCRIPTION:
%   This function computes the optimal pairing between vectors in matrices A and B
%   that minimizes the sum of Euclidean distances. It uses the Hungarian algorithm
%   to solve the assignment problem. If A and B have different numbers of vectors,
%   dummy vectors are added to the smaller matrix to balance the problem.
%
% REQUIRES:
%   MATLAB's matchpairs function (Optimization Toolbox)
%
% EXAMPLE:
%   A = rand(3,5); % 5 vectors in 3D space
%   B = rand(3,5); % 5 vectors in 3D space
%   pairs = optimal_pairing(A, B);
%   % pairs will contain the optimal matching between vectors in A and B

% Initialize variables to store matrix dimensions
m = [0,0];
[N1,m(1)] = size(A);  % Get dimensions of matrix A: N1 rows, m(1) columns
[N2,m(2)] = size(B);  % Get dimensions of matrix B: N2 rows, m(2) columns

% Check if the vectors in A and B have the same dimensionality (number of rows)
if N1 ~= N2
    error('The two input matrices should have the same number of rows.');
else
    N = N1;  % Store the common dimensionality
end

% Check if the two matrices have the same number of vectors (columns)
if m(1) == m(2)
    M = m(1);  % If equal, store the common number of vectors
else
    % If the number of vectors in A and B are different, add dummy
    % vectors (columns) to the smaller matrix to balance the assignment problem
    [M,idx] = max(m);  % M = maximum number of columns, idx = index of the larger matrix
    
    % Create a dummy vector with very large values to ensure it won't be part of
    % the optimal solution unless necessary
    dummy_vect = 1e5+max(max(abs(A),[],'all'),max(abs(B),[],'all'))^2 + zeros(N,1);
    n_dummy = m(idx)-m(setdiff(1:2,idx));  % Calculate how many dummy vectors needed
    
    % Add dummy vectors to the smaller matrix
    if idx == 2  % B has more columns, so add dummy vectors to A
        A = [A, repmat(dummy_vect,1,n_dummy)];
    elseif idx == 1  % A has more columns, so add dummy vectors to B
        B = [B, repmat(dummy_vect,1,n_dummy)];
    end
end

% Compute the pairwise Euclidean distance matrix
% D(i,j) represents the distance between vector A(:,i) and B(:,j)
D = zeros(M, M);
for i = 1:M
    for j = 1:M
        D(i, j) = norm(A(:,i) - B(:,j));  % Calculate Euclidean distance
    end
end

% Solve the assignment problem using the Hungarian algorithm
% The algorithm finds the optimal pairing that minimizes the total distance
% The cost matrix D represents the distances between all possible pairs
% The 1000 parameter is a large value used to replace Inf in the cost matrix
% 'min' indicates we are minimizing the total cost
assignment = matchpairs(D, 1000, 'min');

% Note: assignment is an Mx2 matrix where each row [i,j] means
% vector A(:,i) is paired with vector B(:,j) in the optimal assignment
end