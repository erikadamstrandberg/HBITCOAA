%% Constants
% Fundamental
hbar = 6.62607004e-34/(2*pi);
q_e = 1.60217662e-19;
m_e = 9.109e-31;
c = 2.9979e8;
alpha = 1/137;

% Unit conversion
eV = 1.60217662e-19;
millibarn_per_sr = 1e31/(2*pi);

% Functions
gamma = @(E) E/(m_e*c^2);
v = @(E) c*sqrt(1 - 1./gamma(E).^2);
beta = @(E) v(E)/c;

% Experiment
Z_Ca40 = 20;
Z_e = 1;
E_250MeV = 250e6*eV;
theta = load("data/theta.tsv");
theta_rad = theta*pi/180;
cross_section = load("data/cross_section.tsv");
meas_error = load("data/meas_error.tsv");

%% Relativitic electron?
% For a classical electron gamma approx 1
gamma_experiment = gamma(E_250MeV);

% Assuming all of the energy goes to the E_kin of the electron
v_e_experiment = v(E_250MeV)/c;

%% Proton density and fitting vector
% Proton density that we want to optimize
% X(1) = rho_0, X(2) = a, X(3) = b
rho_ch = @(X, r) X(1)./(1 + exp((r-X(2))./X(3)));

% Create initial guesses for the fitting!
% rho_0 highest proton density
% a middle of fall off
% b sharpness of the falls off
rho_0 = q_e*0.08e45;
a_0 = 4e-15;  
b_0 = 0.1e-15;
X_0 = [rho_0, a_0, b_0];

% Plot of the proton desnsity with the initial guesses
figure(1)
clf; hold on;
r = linspace(0,8,1000)*1e-15;
plot(r, rho_ch(X_0, r))

%% Theoretical cross-section
% Relativistic kinetic energy for a electron
p_e = @(E) sqrt((E/c).^2 - m_e^2*c^2); 

% Difference of p_e before and after collision with Ca40
% Assumes the heavy core to be stationary and not recoil from collision 
q_diff = @(E, theta) sqrt(4*p_e(E).^2.*sin(theta/2).^2);

% Form factor integral broken into two pieces
form_const      = @(E, theta) 4*pi*hbar./(Z_Ca40*q_diff(E,theta)*q_e);
integrand       = @(X, r, E, theta) r.*rho_ch(X, r).*sin(q_diff(E, theta).*r/hbar);
form_integral   = @(X, E, theta) integral(@(r) integrand(X, r, E, theta), 0, 100e-10,'ArrayValued',true);
formfactor      = @(X, E, theta) form_integral(X, E, theta).*form_const(E_250MeV, theta);

% The Rutherford and Mott cross-sections
cross_ruth      = @(E, theta) (Z_e*Z_Ca40*alpha*hbar*c)^2./(4*beta(E)^4*E^2*sin(theta/2).^4);
cross_mott      = @(E, theta) cross_ruth(E, theta).*(1 - beta(E)^2*sin(theta/2).^2);

% Together they make the theoretical cross-section to be fitted against the
% measurement!
cross_theo      = @(E, theta, X) cross_mott(E_250MeV, theta_rad).*abs(formfactor(X_0,E_250MeV,theta)).^2;

figure(2)
clf; hold on;
plot(theta_rad, log10(cross_theo(E_250MeV,theta_rad,X_0)*millibarn_per_sr),'black')
plot(theta_rad, log10(cross_section),'redx')

%%
% Optimization
eps = 1e-6;
integrand       = @(r, X) rho_ch(X, r).*r.^2;
total_charge    = @(X, E) (4*pi/p_e(E))*integral(@(r) integrand(X, r), 0, 1e-10);
constraint      = @(X, E) ((Z_Ca40 - total_charge(X, E))/eps).^2;
Xi2             = @(X, theta) sum(((cross_theo(E_250MeV, theta, X) - cross_section)./meas_error).^2,2);

f = @(X, theta) Xi2(X, theta) + constraint(X, E_250MeV);
X_fit = lsqcurvefit(f, X_0, theta_rad, hannas_very_good_idea);








