function X = kr(A, B)
%======================================================================
% Function: kr
% Purpose : Compute the Khatri–Rao product (column-wise Kronecker product)
%
% Inputs  :
%   A - First input matrix of size (I × K)
%   B - Second input matrix of size (J × K)
%
% Output  :
%   X - Khatri–Rao product of A and B, with size (I*J × K)
%
% Description:
%   The Khatri–Rao product is a column-wise Kronecker product. For each
%   column index k, the function computes:
%
%       X(:,k) = kron(A(:,k), B(:,k))
%
%   This means that the k-th column of the output X is obtained by taking
%   the Kronecker product of the k-th columns of A and B.
%
%   Mathematically:
%       X = [ kron(A(:,1),B(:,1)), kron(A(:,2),B(:,2)), ..., kron(A(:,K),B(:,K)) ]
%
%   This operation is widely used in tensor decomposition (e.g., PARAFAC)
%   and multi-way signal processing.
%======================================================================

[~, K] = size(B);        % Get the number of columns (assumed same for A and B)
for k = 1:K              % Loop over each column index
    X(:, k) = kron(A(:, k), B(:, k));  % Compute column-wise Kronecker product
end
%======================================================================
end
