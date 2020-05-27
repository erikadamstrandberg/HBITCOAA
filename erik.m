%% Constants
global hbar q_e m_e c alpha
% Fundamental
hbar = 6.62607004e-34/(2*pi);
q_e = 1.60217662e-19;
m_e = 9.109e-31;
c = 2.9979e8;
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
data = [theta_rad, cross_section, meas_error, Z_Ca40, Z_e, E_250MeV];

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
rho_0 = q_e*0.07e45;
a_0 = 3.44e-15;  
b_0 = 0.2e-15;
X_0 = [rho_0, a_0, b_0];

% Plot of the proton desnsity with the initial guesses
figure(1)
clf; hold on;
r = linspace(0,100,1000)*1e-15;
plot(r, rho_ch(X_0, r))

%% Theoretical cross-section
% Relativistic kinetic energy for a electron
p_e             = @(E) sqrt((E/c).^2 - m_e^2*c^2); 

% Difference of p_e before and after collision with Ca40
% Assumes the heavy core to be stationary and not recoil from collision 
q_diff          = @(E, theta) sqrt(4*p_e(E).^2.*sin(theta/2).^2);

upper_limit = 18e-15;
% Form factor integral broken into two pieces
form_const      = @(E, theta) 4*pi*hbar./(Z_Ca40*q_diff(E,theta)*q_e);
integrand       = @(X, r, E, theta) r.*rho_ch(X, r).*sin(q_diff(E, theta).*r/hbar);
form_integral   = @(X, E, theta) integral(@(r) integrand(X, r, E, theta), 0, upper_limit,'ArrayValued',true);
formfactor      = @(X, E, theta) form_integral(X, E, theta).*form_const(E_250MeV, theta);

% The Rutherford and Mott cross-sections
cross_ruth      = @(E, theta) (Z_e*Z_Ca40*alpha*hbar*c)^2./(4*beta(E)^4*E^2*sin(theta/2).^4);
cross_mott      = @(E, theta) cross_ruth(E, theta).*(1 - beta(E)^2*sin(theta/2).^2);

% Together they make the theoretical cross-section to be fitted against the
% measurement!
cross_theo      = @(E, theta, X) cross_mott(E, theta).*abs(formfactor(X_0,E,theta)).^2;

figure(2)
clf; hold on;
plot(theta, log10(cross_theo(E_250MeV, theta_rad, X_0)*millibarn_per_sr),'black')
plot(theta, log10(cross_section),'redx')


%% Theoretical cross-section from box distribution

box_integrand  = @(X, r, E, theta) r.*X(1).*sin(q_diff(E, theta).*r/hbar);
box_integral   = @(X, E, theta) integral(@(r) box_integrand(X, r, E, theta), 0, X(2),'ArrayValued',true);
box_formfactor = @(X, E, theta) box_integral(X, E, theta).*form_const(E_250MeV, theta);
cross_box      = @(E, theta, X) cross_mott(E, theta).*abs(box_formfactor(X_0,E,theta)).^2;

theta_plot = linspace(0,180,1000)';
theta_plot_rad = theta_plot*pi/180;

figure(3)
clf; hold on;
plot(theta_plot, log10(cross_theo(E_250MeV,theta_plot_rad, X_0)),'black')
plot(theta_plot, log10(cross_box(E_250MeV,theta_plot_rad, X_0)),'red')
legend("rho saxton","rho box")

%% Theoretical cross-section with pertubated energy

theta_plot = linspace(35,135,1000)';
theta_plot_rad = theta_plot*pi/180;

pertubation = 10e6*eV;
delta_eV = 5e6*eV;
energies = (E_250MeV - pertubation):delta_eV:(E_250MeV + pertubation);

figure(4)
clf; hold on;
for i = 1:length(energies)
    plot(theta_plot, log10(cross_theo(energies(i),theta_plot_rad, X_0)*millibarn_per_sr))
end


%%
% Plotting the formfactor integrand to that it is okey with a upper limit
integrand_plot = integrand(X_0, r, E_250MeV, theta_rad);
x_plot = [upper_limit, upper_limit];
y_plot = [min(min(integrand_plot)), max(max(integrand_plot))];
max_value = max(integrand(X_0, upper_limit, E_250MeV, theta_rad));

figure(5)
clf; hold on;
plot(r, integrand_plot)
plot(x_plot,y_plot,'--')
axis([0,20e-15,y_plot(1),y_plot(2)])
%%
% Optimization
% eps = 1e-6;
% integrand       = @(r, X) rho_ch(X, r).*r.^2;
% total_charge    = @(X, E) (4*pi/p_e(E))*integral(@(r) integrand(X, r), 0, 1e-10);
% constraint      = @(X, E) ((Z_Ca40 - total_charge(X, E))/eps).^2;
% Xi2             = @(X, theta) sum(((cross_theo(E_250MeV, theta, X) - cross_section)./meas_error),2);
% 
% f = @(X, theta) Xi2(X, theta) + constraint(X, E_250MeV);
% X_fit = lsqcurvefit(f, X_0, theta_rad, hannas_very_good_idea);


