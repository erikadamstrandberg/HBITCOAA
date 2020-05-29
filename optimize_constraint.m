function F = optimize_constraint(X_0, rho_ch, exp_values, q_e)

Z_Ca40 = exp_values(1);
E_250MeV = exp_values(3);

integrand       = @(X, r) rho_ch(X, r).*r.^2;
total_charge    = @(X) (4*pi/q_e)*integral(@(r) integrand(X, r), 0, 100e-15);
constraint      = @(X) ((Z_Ca40 - total_charge(X))).^2;

F = constraint(X_0);