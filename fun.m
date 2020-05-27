function F = fun(X_0, rho_ch, p_e, cross_theo, data, exp_values, millibarn_per_sr)

eps = 1e-6;
integrand       = @(X, r) rho_ch(X, r).*r.^2;
total_charge    = @(X, E) (4*pi/p_e(E))*integral(@(r) integrand(X, r), 0, 1e-10);
constraint      = @(X) ((exp_values(1) - total_charge(X, exp_values(3)))/eps).^2;
Xi2             = @(X) ((cross_theo(exp_values(3), data(1), X)*millibarn_per_sr - data(2))./data(3)).^2;

F = Xi2(X_0) + ones(length(data(1)), 1)*constraint(X_0);