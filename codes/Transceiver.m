function [trans_X, trans_X_inv] = Transceiver(M, T)
%======================================================================
% Function: Transceiver
% Purpose : Generate the transmitted pilot signal matrix (X) used at the BS
% Inputs  :
%   M - Number of transmitting antennas (or users)
%   T - Number of time slots (pilot symbols)
% Outputs :
%   trans_X     - The generated M×T pilot signal matrix
%   trans_X_inv - The pseudo-inverse (or matched filter) of trans_X
%
% Description:
%   This function generates a structured transmit pilot matrix based on
%   the Discrete Fourier Transform (DFT). The generated pilot signals are
%   orthogonal and have constant modulus, ensuring low correlation between
%   pilot sequences.
%
%   Depending on whether M < T or M ≥ T, the DFT matrix is truncated or
%   resized appropriately. The inverse (trans_X_inv) is computed based on
%   the corresponding matrix dimensions.
%
%   The commented section below provides an alternative approach using
%   random SVD-based orthogonal matrices.
%======================================================================

% DFT-based orthogonal pilot generation (active) ---
if M < T
    % Case 1: Fewer antennas than time slots
    TEMP = dftmtx(T);                         % Generate T×T DFT matrix
    trans_X = 1/sqrt(T) * TEMP(1:M, 1:T);     % Select first M rows
    trans_X_inv = trans_X' * inv(trans_X * trans_X');  % Compute pseudo-inverse
else
    % Case 2: More antennas than time slots
    TEMP = dftmtx(M);                         % Generate M×M DFT matrix
    trans_X = 1/sqrt(M) * TEMP(1:M, 1:T);     % Select first T columns
    trans_X_inv = inv(trans_X' * trans_X) * trans_X';  % Compute pseudo-inverse
end
%======================================================================
end
