function [holo] = holoNorm(holoRaw, bkgrnd)

% 11-Oct-2017 JF
%
% holoNorm: function that normalises a hologram by a background.
%
% Input:    raw hologram 
%           background
%           *Must be same size and double precision*
%
% Outputs:  normalised hologram
%           row dimension of new array
%           col dimension of new array
%
% Normalisation proceeds:
% 
%   normalisedHolo = (rawHolo - background)/(2*sqrt(background))
%
% Reference: Lee & Grier, Optics Express 15 (2007)
% ===================================

%% Subtract background
sub = holoRaw-bkgrnd;

%% Divide
holo = sub./(2*sqrt(bkgrnd));

