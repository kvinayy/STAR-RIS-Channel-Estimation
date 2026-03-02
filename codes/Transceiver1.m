function [trans_X, trans_X_inv] = Transceiver1(M1, T)
%======================================================================
% Function: Transceiver1
% Purpose : Generate the transmitted pilot signal matrix (X) for the 
%           transmission-side users in the STAR-RIS system.
%
% Inputs  :
%   M1 - Number of transmission-side users (or transmit antennas)
%   T  - Number of time slots (pilot symbols)
%
% Outputs :
%   trans_X     - The generated M1×T pilot signal matrix
%   trans_X_inv - The pseudo-inverse (or matched filter) of trans_X
%
% Description:
%   This function generates the pilot transmission matrix used for the
%   transmission (T) users of the STAR-RIS system. It creates structured
%   orthogonal pilot sequences using the Discrete Fourier Transform (DFT)
%   matrix, which provides good orthogonality and constant modulus.
%
%   Depending on whether M1 < T or M1 ≥ T, the DFT matrix is generated and
%   truncated accordingly to match the required size. The corresponding
%   pseudo-inverse is then calculated to allow signal reconstruction.
%
%   The commented section below shows an alternative random SVD-based
%   orthogonal pilot generation method.
%======================================================================

% DFT-based orthogonal pilot generation (active) ---
if M1 < T
    % Case 1: Fewer users (M1) than time slots (T)
    TEMP = dftmtx(T);                          % Generate T×T DFT matrix
    trans_X = 1/sqrt(T) * TEMP(1:M1, 1:T);     % Select first M1 rows
    trans_X_inv = trans_X' * inv(trans_X * trans_X');  % Compute pseudo-inverse
else
    % Case 2: More users (M1) than time slots (T)
    TEMP = dftmtx(M1);                         % Generate M1×M1 DFT matrix
    trans_X = 1/sqrt(M1) * TEMP(1:M1, 1:T);    % Select first T columns
    trans_X_inv = inv(trans_X' * trans_X) * trans_X';  % Compute pseudo-inverse
end
%======================================================================
end
