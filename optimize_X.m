function F = optimize_X(X_0, rho_ch, cross_theo, data, exp_values, millibarn_per_sr, q_e)

theta_rad = data(:,1);
cross_section = data(:,2);
meas_error = data(:,3);

Z_Ca40 = exp_values(1);
E_250MeV = exp_values(3);

eps = 1e-6;
integrand       = @(X, r) rho_ch(X, r).*r.^2;
total_charge    = @(X) (4*pi/q_e)*integral(@(r) integrand(X, r), 0, 100e-15);
constraint      = @(X) ((Z_Ca40 - total_charge(X))/eps).^2;
Xi2             = @(X) ((cross_theo(X, E_250MeV, theta_rad)*millibarn_per_sr - cross_section)./meas_error).^2;

F = [Xi2(X_0); constraint(X_0)];