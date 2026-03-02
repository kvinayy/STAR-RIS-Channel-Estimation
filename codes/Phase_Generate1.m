function [Theta] = Phase_Generate1(P, N)
%======================================================================
% Function: Phase_Generate1
% Purpose : Generate the pilot phase matrix (Θ) for the STAR-RIS system
% Inputs  : 
%   P - Number of rows (number of pilot configurations or phase patterns)
%   N - Number of RIS elements
% Output  : 
%   Theta - Generated P×N pilot phase matrix
%
% Description:
%   This function generates the RIS phase-shift matrix Θ used for the 
%   transmission side of the STAR-RIS during training.
%
%   The phase matrix can be generated either randomly or using a 
%   DFT-based structured approach.
%
%   - The commented lines below show alternative generation methods:
%       1) Random phase matrix using complex exponential
%       2) Using the first N columns of a unitary matrix obtained from SVD
%
%   - The active implementation uses the Discrete Fourier Transform (DFT)
%     matrix to form structured and orthogonal phase patterns.
%======================================================================

% DFT-based phase generation (currently used) ---
TEMP = dftmtx(N);                 % Generate the N×N DFT matrix
% Theta = 1/sqrt(N) * TEMP(1:P, 1:N);   % (Normalized DFT-based version)
Theta = sqrt(0.5) * TEMP(1:P, 1:N);    % (DFT-based version with scaling factor)
%======================================================================
end
