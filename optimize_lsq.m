function F = optimize_lsq(X_0, theta, rho_ch, cross_theo, data, exp_values, millibarn_per_sr, q_e)

meas_error = data(:,3);
Z_Ca40 = exp_values(1);
E_250MeV = exp_values(3);

eps = 1e-6;
integrand       = @(X, r) rho_ch(X, r).*r.^2;
total_charge    = @(X) (4*pi/q_e)*integral(@(r) integrand(X, r), 0, 100e-15);
constraint      = @(X) (Z_Ca40 - total_charge(X))/eps;
Xi2             = @(X, theta, meas_error) cross_theo(X, E_250MeV, theta)*millibarn_per_sr./meas_error;

F = zeros(length(theta),1);
for i = 1:length(F)
    if i == 32
        F(i) = constraint(X_0);
    else
        F(i) = Xi2(X_0, theta(i), meas_error(i));
    end
end