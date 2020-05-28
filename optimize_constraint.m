function F = optimize_constraint(X_0, rho_ch, exp_values, q_e, eps)

Z_Ca40 = exp_values(1);
E_250MeV = exp_values(3);

integrand       = @(X, r) rho_ch(X, r).*r.^2;
total_charge    = @(X, E) (4*pi/q_e)*integral(@(r) integrand(X, r), 0, 20e-15);
constraint      = @(X) ((Z_Ca40 - total_charge(X, E_250MeV))/eps).^2;

F = constraint(X_0);