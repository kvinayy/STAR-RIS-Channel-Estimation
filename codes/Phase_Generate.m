function [Phi] = Phase_Generate(P, N)
%======================================================================
% Function: Phase_Generate
% Purpose : Generate the pilot phase matrix (Φ) for the STAR-RIS system
% Inputs  : 
%   P - Number of rows (number of pilot configurations or phase patterns)
%   N - Number of RIS elements
% Output  : 
%   Phi - Generated P×N pilot phase matrix
%
% Description:
%   This function generates the RIS phase-shift matrix Φ used for training.
%   It can be generated either randomly or using a DFT-based approach.
%
%   - The commented lines below show alternative methods:
%       1) Random phase matrix using exp(j·2π·rand)
%       2) Using the first N columns of a unitary matrix from SVD
%
%   - The active implementation uses the Discrete Fourier Transform (DFT)
%     matrix to construct structured phase patterns.
%======================================================================

% DFT-based phase generation ---
TEMP = dftmtx(N);               % Generate the N×N DFT matrix
% Phi = 1/sqrt(N) * TEMP(1:P, 1:N);  % (Normalized DFT-based version)
Phi = sqrt(0.5) * TEMP(1:P, 1:N);   % (DFT-based version with scaling factor)
%======================================================================
end
