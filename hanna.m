%% Constants
% Fundamental
hbar = 6.62607004e-34/(2*pi);
q_e = 1.60217662e-19;
m_e = 9.109e-31;
eV = 1.60217662e-19;
c = 2.9979e8;
alpha = 1/137;

% Experiment
Z_Ca40 = 20;
Z_e = 1;
E_experiment = 250e6*eV;         % (250 MeV)
theta = load("data/theta.tsv");
cross_section = load("data/cross_section.tsv");
meas_error = load("data/meas_error.tsv");


%%
% functions
rho_ch = @(X, r) X(1)./(1 + exp((r-X(2))./X(3)));
p = @(E, m) sqrt((E/c).^2 - (m*c)^2);
gamma = @(E, m) E/(m*c^2);
v = @(E, m) c*sqrt(1 - 1./gamma(E, m).^2);
beta = v(E_experiment, m_e)/c;


p_e = p(E_experiment, m_e);
v_e = v(E_experiment, m_e);


%%
cross_ruth = @(E, theta) (Z_e*Z_Ca40*alpha*hbar*c)^2./(4*beta^4*E^2*sin(theta/2).^4);
cross_mott = @(E, theta) cross_ruth(E, theta).*(1 - beta^2*sin(theta/2).^2);

cross_theo = @(E, theta, X) cross_mott(E, theta)*(X(1)+X(2)+X(3)); % for test only


theta_vect = linspace(0, pi, 100);
cross_section_theo = cross_mott(E_experiment, theta_vect);


clf
hold on 
plot(theta_vect, log(cross_section_theo))


%%
% Constraint
integrand = @(X, r) rho_ch(X, r).*r.^2;
total_charge = @(X) (4*pi/p_e)*integral(@(r) integrand(X, r), 0, 1e-10);

X_0 = [q_e*0.073463e45, 4e-15, 0.1e-15];
total_charge(X_0)

%%
eps = 1e-6;
constraint = @(X) ((Z_Ca40 - total_charge(X))/eps).^2;

R = linspace(0, 5e-15, 100); 


clf
plot(R, integrand(X_0, R));

constraint(X_0)

%%
Xi2 = @(X, theta) ((cross_theo(E_experiment, theta, X) - cross_section)./meas_error).^2;

Xi2(X_0,theta_rad)

% f = @(X, theta) Xi2(X, theta) + constraint(X);
% 
% X_fit = lsqcurvefit(f, X_0, theta.', cross_section.');

%%

%% HANNAS OLD

%% %%%%%%%%%%%  Optimization %%%%%%%%%%%%
% optimize the constraint
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off');
X_con = fsolve(@(X) optimize_constraint(X, rho_ch, exp_values, q_e, 1), X_0, options);

% Check charge
integrand       = @(X, r) rho_ch(X, r).*r.^2;
total_charge    = @(X, E) (4*pi/q_e)*integral(@(r) integrand(X, r), 0, 20e-15);

charge = total_charge(X_con, E_250MeV);
fprintf("Total charge: " + charge + "\n")

% optimize cross-section
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off');
X_star = fsolve(@(X) optimize_X(X, rho_ch, cross_theo, data, exp_values, millibarn_per_sr, q_e), X_con, options);


%%
options_less_eps = optimoptions('fsolve','Algorithm','levenberg-marquardt', 'Display', 'iter', 'TolX', 1e-12);

X_con = fsolve(@(X) optimize_constraint(X, rho_ch, exp_values, q_e, 1e-6), X_iter, options_less_eps);

%%
X_star = fsolve(@(X) optimize_X(X, rho_ch, cross_theo, data, exp_values, millibarn_per_sr, q_e), X_con, options_less_eps);




