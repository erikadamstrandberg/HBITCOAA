
%% Constants
% Fundamental
hbar = 6.62607004e-34/(2*pi);
q_e = 1.60217662e-19;
m_e = 9.109e-31;
c = 2.9979e8;
alpha = 1/137;

% Unit conversion
eV = 1.60217662e-19;
millibarn = 1e31;

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
% Rho_0
% a is in the middle of the drop
% b is how sharp the drop falls off
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
p_e = @(E) sqrt(E/c^2 - m_e^2*c^2);

% Difference of p_e before and after collision with Ca40
% Assumes the heavy core to be stationary and not recoil from collision 
q_diff = @(E, theta) 4*p_e(E).^2.*sin(theta/2).^2;

% Form factor integral broken into two pieces
form_const = @(E, theta) 4*pi*hbar./(Z_Ca40*q_diff(E,theta));
integrand = @(X, r, E, theta) r.*rho_ch(X, r).*sin(sqrt(q_diff(E, theta).*r/hbar));
form_integral = @(X, E, theta) integral(@(r) integrand(X, r, E, theta), 0, 1000e-15);

n = length(theta_rad);
formfactor = zeros(n, 1);
for i = 1:n
    formfactor(i) = form_integral(X_0, E_250MeV, theta_rad(i))*form_const(E_250MeV, theta_rad(i));
end

cross_ruth = @(E, theta) (Z_e*Z_Ca40*alpha*hbar*c)^2./(4*beta(E)^4*E^2*sin(theta/2).^4);
cross_mott = @(E, theta) cross_ruth(E, theta).*(1 - beta(E)^2*sin(theta/2).^2);
cross_theo = cross_mott(E_250MeV, theta_rad).*abs(formfactor).^2;

figure(2)
clf; hold on;
plot(theta_rad, log10(cross_theo*millibarn))






