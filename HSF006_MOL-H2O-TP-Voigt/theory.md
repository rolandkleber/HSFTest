Main Theory pdf https://doi.org/10.5281/zenodo.16953782
% theory.m — MOL-01 · H2O Linewidth T/P Sweep (Voigt)
% Minimal theory & quick calculators (MATLAB/Octave)
%
% Parameters inferred from last run (summary.json if present):
%   nu0   = 7181.155657 cm^-1
%   gamma_ref (HWHM at Tref,Pref) = 0.099700 cm^-1
%   n_exp (temperature exponent)   = 0.760
%   Tref  = 296.0 K, Pref = 1.0 atm
%
% Example:
%   T = [250 300 350 400 500 700];
%   P = [0.2 0.5 1.0];
%   [FWHM_G, gammaL] = h2o_widths(T, P, 7181.155657, 0.099700, 0.760, 296.0, 1.0);
%   slope375 = expected_dgamma_dP(375, 0.099700, 0.760, 296.0);
%
% Notes:
%   - FWHM_G is pure Doppler (no instrument).
%   - gammaL is air-broadened collisional HWHM.

function [FWHM_G_cm1, gammaL_cm1] = h2o_widths(T_K, P_atm, nu0_cm1, gamma_ref_cm1, n_exp, Tref_K, Pref_atm)
  % Physical constants
  kB = 1.380649e-23;   % J/K
  c  = 299792458.0;    % m/s
  amu= 1.66053906660e-27; % kg
  mH2O = 18.01528 * amu;  % kg

  T_K = T_K(:)';
  P_atm = P_atm(:)';

  % Doppler FWHM (cm^-1)
  FWHM_G_cm1 = nu0_cm1 .* sqrt(8*kB*log(2).*T_K./(mH2O*c^2));

  % Collisional HWHM (cm^-1)
  gammaL_cm1 = gamma_ref_cm1 .* (Tref_K ./ T_K).^n_exp .* (P_atm./Pref_atm);
end

function slope = expected_dgamma_dP(T_K, gamma_ref_cm1, n_exp, Tref_K)
  % Expected linear slope d(gamma)/dP at temperature T_K (cm^-1 per atm)
  slope = gamma_ref_cm1 * (Tref_K / T_K)^n_exp;
end
