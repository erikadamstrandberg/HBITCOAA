% ----------------------------------------------------------------------
% Creates the function to be optimized with lsqcurvefit
% ----------------------------------------------------------------------
function F = optimize_lsq(X_0, xdata, rho_ch, cross_theo, data, exp_values, eps_charge, millibarn, q_e)
% Arguments:
%
%   X_0                     = initial guesses         [rho_0, a_0, b_0]
%   xdata                   = x-data for fitting      [theta_rad; 0]
%   rho_ch(X, r)            = calculates charge distribution
%   cross_theo(X, E, theta) = calculates theoretical cross-section
%   data                    = data from measurement
%   exp_values              = values used in experiment
%   eps_charge              = Weight for the conservation of charge constraint

% Unpack measurement and experiment values
meas_error = data(:,3);
theta_rad = xdata(1:length(meas_error));
E_250MeV = exp_values(3);

integrand       = @(X, r) rho_ch(X, r).*r.^2;
total_charge    = @(X) (4*pi/q_e)*integral(@(r) integrand(X, r), 0, 100e-15);
constraint      = @(X) total_charge(X)/eps_charge;
Xi2             = @(X, theta, meas_error) cross_theo(X, E_250MeV, theta)./(meas_error*millibarn);

% Note that the last value is the conservation of charge and that it is
% independet of theta or the xdata for the fit.
F = [Xi2(X_0, theta_rad, meas_error); constraint(X_0)];