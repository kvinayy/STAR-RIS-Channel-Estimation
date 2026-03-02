function [a] = steering_vector(array,theta, phi)
%======================================================================
% Function: steering_vector
% Purpose : Compute the array steering vector for a given direction
%
% Inputs  :
%   array - Matrix containing 3D positions of antenna elements
%           (each row corresponds to one antenna element)
%   theta - Azimuth angle of arrival (radians)
%   phi   - Elevation angle of arrival (radians)
%
% Output :
%   a     - Complex steering vector representing phase shifts
%           across the antenna array
%
% Description:
%   This function models the phase response of each antenna element
%   when a plane wave arrives from direction (theta, phi).
%   The phase shift depends on the dot product of antenna position
%   and the wave vector.
%======================================================================

    % Compute steering vector using spatial phase shifts
a = exp(-1i*array*K(theta,phi));
end

function k = K(theta,phi)
%======================================================================
% Function: K
% Purpose : Compute the wave (spatial frequency) vector
%
% Inputs  :
%   theta - Azimuth angle (radians)
%   phi   - Elevation angle (radians)
%
% Output :
%   k     - 3x1 wave vector corresponding to the signal direction
%
% Description:
%   The wave vector defines the propagation direction of the signal.
%   The factor pi assumes half-wavelength antenna spacing.
%======================================================================

    % Direction-dependent wave vector
k = pi * [cos(theta).*cos(phi), sin(theta).*cos(phi), sin(phi)]';
end

