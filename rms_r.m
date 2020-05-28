function radius = rms_r(X_star, rho_ch, exp_values, q_e)

integrand = @(X, r) r.^4.*rho_ch(X, r);
radius_squared = (4*pi/(exp_values(1)*q_e))*integral(@(r) integrand(X_star, r), 0, 100e-10);

radius = sqrt(radius_squared);