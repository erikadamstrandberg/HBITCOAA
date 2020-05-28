%% Constants
global hbar q_e m_e c alpha eV millibarn_per_sr
% Fundamental
hbar  = 6.62607004e-34/(2*pi);
q_e   = 1.60217662e-19;
m_e   = 9.109e-31;
c     = 2.9979e8;
alpha = 1/137;

% Unit conversion
eV = 1.60217662e-19;
millibarn_per_sr = 1e31*2*pi;

% Functions
gamma = @(E) E/(m_e*c^2);
v = @(E) c*sqrt(1 - 1./gamma(E).^2);
beta = @(E) v(E)/c;

% Measurement data
theta = load("data/theta.tsv");
theta_rad = theta*pi/180;
cross_section = load("data/cross_section.tsv");
meas_error = load("data/meas_error.tsv");

% Experiment
Z_Ca40 = 20;
Z_e = 1;
E_250MeV = 250e6*eV;

% Used for passing data to functions
data = [theta_rad, cross_section, meas_error];
exp_values = [Z_Ca40, Z_e, E_250MeV];

%% Proton density and fitting vector
scale_rho = (1/1.2817)*1e-25;
scale_a = (1/4)*1e15;
scale_b = (1/1.6)*1e15;
scale = [scale_rho, scale_a, scale_b];

% Proton density that we want to optimize
% X(1) = rho_0, X(2) = a, X(3) = b
rho_ch = @(X, r) (1/scale(1))*X(1)./(1 + exp((r-(1/scale(2))*X(2))./(X(3)*(1/scale(3)))));
%rho_ch = @(X, r) (1/scale_factor)*X(1)./(1 + exp((r-X(2))./X(3)));

% Create initial guesses for the fitting!
% rho_0 highest proton density
% a middle of fall off
% b sharpness of the falls off
rho_0 = q_e*0.08e45*scale(1); 
a_0 = 4e-15*scale(2);  
b_0 = 1.6e-15*scale(3);
X_0 = [1, 0.9, 0.5]; %[rho_0, a_0, b_0];

% Plot of the proton desnsity with the initial guesses
figure(1)
clf; hold on;
r = linspace(0,8,1000)*1e-15;
plot(r, rho_ch(X_0, r))

% Theoretical cross-section
% Relativistic kinetic energy for a electron
p_e = @(E) sqrt((E/c).^2 - m_e^2*c^2); 

% Difference of p_e before and after collision with Ca40

% Assumes the heavy core to be stationary and not recoil from collision 
q_diff = @(E, theta) sqrt(4*p_e(E).^2.*sin(theta/2).^2);

% Form factor integral broken into two pieces
form_const      = @(E, theta) 4*pi*hbar./(Z_Ca40*q_diff(E,theta)*q_e);
form_integrand  = @(X, r, E, theta) r.*rho_ch(X, r).*sin(q_diff(E, theta).*r/hbar);
form_integral   = @(X, E, theta) integral(@(r) form_integrand(X, r, E, theta), 0, 100e-10,'ArrayValued',true);
formfactor      = @(X, E, theta) form_integral(X, E, theta).*form_const(E_250MeV, theta);

% The Rutherford and Mott cross-sections
cross_ruth      = @(E, theta) (Z_e*Z_Ca40*alpha*hbar*c)^2./(4*beta(E)^4*E^2*sin(theta/2).^4);
cross_mott      = @(E, theta) cross_ruth(E, theta).*(1 - beta(E)^2*sin(theta/2).^2);

% Together they make the theoretical cross-section to be fitted against the
% measurement!
cross_theo      = @(E, theta, X) cross_mott(E_250MeV, theta_rad).*abs(formfactor(X_0,E_250MeV,theta)).^2;

figure(2)
clf; hold on;
plot(theta_rad, log10(cross_theo(E_250MeV,theta_rad,X_0)),'black')
plot(theta_rad, log10(cross_section/millibarn_per_sr),'redx')

%% %%%%%%%%%%%  Optimization %%%%%%%%%%%%
% optimize the constraint
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter');
X_con = fsolve(@(X) optimize_constraint(X, rho_ch, exp_values, q_e, 1), X_0, options);

% Check charge
integrand       = @(X, r) rho_ch(X, r).*r.^2;
total_charge    = @(X, E) (4*pi/q_e)*integral(@(r) integrand(X, r), 0, 20e-15);

charge = total_charge(X_con, E_250MeV);
fprintf("Total charge: " + charge + "\n")

% optimize cross-section
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter');
X_star = fsolve(@(X) optimize_X(X, rho_ch, cross_theo, data, exp_values, millibarn_per_sr, q_e), X_con, options);


%% Iterate between opt. constrain and opt. cross-section
niter = 2;
X_iter = X_0;

options_iter = optimoptions('fsolve','Algorithm','levenberg-marquardt', 'Display', 'iter', 'TolX', 1e-15);

for i=1:niter
    fprintf("Iteration " + i + "\n")
    
    X_con = fsolve(@(X) optimize_constraint(X, rho_ch, exp_values, q_e, 1), X_iter, options_iter);
    X_star = fsolve(@(X) optimize_X(X, rho_ch, cross_theo, data, exp_values, millibarn_per_sr, q_e), X_con, options_iter);
    X_iter = X_star
end

%%
options_less_eps = optimoptions('fsolve','Algorithm','levenberg-marquardt', 'Display', 'iter', 'TolX', 1e-12);

X_con = fsolve(@(X) optimize_constraint(X, rho_ch, exp_values, q_e, 1e-6), X_iter, options_less_eps);

%%
X_star = fsolve(@(X) optimize_X(X, rho_ch, cross_theo, data, exp_values, millibarn_per_sr, q_e), X_con, options_less_eps);

%%
% check optimized theoretic cross-section
figure(3)
clf; hold on;
plot(theta_rad, log10(cross_theo(E_250MeV,theta_rad,X_0)),'black')
plot(theta_rad, log10(cross_section/millibarn_per_sr),'redx')
plot(theta_rad, log10(cross_theo(E_250MeV,theta_rad,X_star)),'blue')