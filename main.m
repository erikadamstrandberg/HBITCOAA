%% Constants
global hbar q_e m_e c alpha eV millibarn_per_sr
% Fundamental
hbar  = 6.62607004e-34/(2*pi);
q_e   = 1.60217662e-19;
m_e   = 9.109e-31;
c     = 299792458;
alpha = 1/137;

% Unit conversion
eV = 1.60217662e-19;
millibarn_per_sr = 1e31*2*pi;

% Functions
gamma = @(E) E/(m_e*c^2);
p_e = @(E) sqrt((E/c).^2 - m_e^2*c^2); 
v = @(E) p_e(E)./(m_e*gamma(E));
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
% Create initial guesses for the fitting!
% rho_0 highest proton density
% a middle of fall off
% b sharpness of the falls off
rho_0 = q_e*0.07*1e45; 
a_0 = 3e-15;  
b_0 = 0.5e-15;

% Scaling for fsolve to work with equally big paramters
X_real = [rho_0, a_0, b_0];
X_scale = [1/rho_0, 1/a_0, 1/b_0];
X_0 = X_real.*X_scale;

% Proton density that we want to optimize
% X(1) = rho_0, X(2) = a, X(3) = b
rho_ch = @(X, r) (X(1)/X_scale(1))./(1 + exp((r-(X(2)/X_scale(2)))./(X(3)/X_scale(3))));

% Plot of the proton desnsity with the initial guesses
figure(1)
clf; hold on;
r = linspace(0,8,1000)*1e-15;
plot(r, rho_ch(X_0, r)/(q_e*1e45))

%%
% Theoretical cross-section
% Difference of p_e before and after collision with Ca40
% Assumes the heavy core to be stationary and not recoil from collision 
q_diff = @(E, theta) sqrt(4*p_e(E).^2*sin(theta/2).^2);

% Form factor integral spilt into two pieces
form_const      = @(E, theta) 4*pi*hbar./(Z_Ca40*q_e*q_diff(E, theta));
form_integrand  = @(X, r, E, theta) r.*rho_ch(X, r).*sin(q_diff(E, theta).*r/hbar);
form_integral   = @(X, E, theta) integral(@(r) form_integrand(X, r, E, theta), 0, 100e-15,'ArrayValued',true);
formfactor      = @(X, E, theta) form_integral(X, E, theta).*form_const(E, theta);

% The Rutherford and Mott cross-sections
cross_ruth      = @(E, theta) (Z_e*Z_Ca40*alpha*hbar*c)^2./(4*beta(E)^4*E^2*sin(theta/2).^4);
cross_mott      = @(E, theta) cross_ruth(E, theta).*(1 - beta(E)^2*sin(theta/2).^2);

% Together they make the theoretical cross-section to be fitted against the
% measurement!
cross_theo      = @(X, E, theta) cross_mott(E, theta).*abs(formfactor(X, E, theta)).^2;

figure(2)
clf; hold on;
plot(theta, log10(cross_theo(X_0, E_250MeV,theta_rad)*millibarn_per_sr),'black')
plot(theta, log10(cross_section),'redx')

%% Iterate between opt. constrain and opt. cross-section

options_constraint = optimoptions('fsolve','Algorithm',           'levenberg-marquardt',...
                                           'Display',             'off',...
                                           'StepTolerance',       1e-15,...
                                           'FunctionTolerance',   1e-10);
X_con = fsolve(@(X) optimize_constraint(X, rho_ch, exp_values, q_e), X_0, options_constraint); 

theta_rad_fit = [theta_rad;0];
cross_section_fit = [cross_section./meas_error;0];

X_star = lsqcurvefit(@(X,theta) optimize_lsq(X, theta, rho_ch, cross_theo, data, exp_values, millibarn_per_sr, q_e),...
                                X_con,...
                                theta_rad_fit,...
                                cross_section_fit);
                           

%%
% check optimized theoretic cross-section
figure(3)
clf; hold on;
plot(theta, log10(cross_theo(X_0, E_250MeV,theta_rad)*millibarn_per_sr),'black')
plot(theta, log10(cross_section),'redx')
plot(theta, log10(cross_theo(X_star, E_250MeV,theta_rad)*millibarn_per_sr),'blue')

radius = rms_r(X_star, rho_ch, exp_values, q_e)
integrand       = @(X, r) rho_ch(X, r).*r.^2;
total_charge    = @(X) (4*pi/q_e)*integral(@(r) integrand(X, r), 0, 100e-15);
total_charge(X_star)

