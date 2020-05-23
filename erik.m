%% Constants
% Fundamental
hbar = 6.62607004e-34/(2*pi);
q_e = 1.60217662e-19;
m_e = 9.109e-31;
eV = 1.60217662e-19;
c = 2.9979e8;

% Functions
v = @(MeV) c*sqrt(1 - 1./gamma(MeV).^2);
E_e = @(MeV) MeV*1e6*eV;

% Experiment
Z_Ca40 = 20;
E_250MeV = E_e(250);
theta = load("data/theta.tsv");
cross_section = load("data/cross_section.tsv");

%% Relativitic electron?
% For a classical electron gamma approx 1
gamma = @(MeV) E_e(MeV)/(m_e*c^2);
gamma_experiment = gamma(250);

% Assuming all of the energy goes to the E_kin of the electron
v_e_experiment = v(250)/c;

%% Theoretical cross-section
% Proton density that we want to optimize
% X(1) = rho_0, X(2) = a, X(3) = b
rho_ch = @(X,r) X(1)./(1 + exp((r-X(2))./X(3)));

% Relativistic kinetic energy for a electron
p_e = @(E_e) sqrt(E_e/c^2 - m_e^2*c^2);

% Difference of p_e before and after collision with Ca40
% Assumes the heavy core to be stationary and not recoil from collision 
q_diff = @(E_e,theta) 4*p_e(E_e).^2.*sin(theta/2).^2;

% Form factor integral broken into two pieces.
form_const = @(E_e,theta) 4*pi*hbar./(Z_Ca40*q_diff(E_e,theta));
form_integrand = @(X,r,E_e,theta) r.*rho_ch(X,r).*sin(sqrt(q_diff(E_e,theta).*r/hbar));
form_integral = @(X,E_e,theta) integral(@(r)form_integrand(X,r,E_e,theta),0,inf);

n = length(theta);
formfactor = zeros(n,1);
for i = 1:n
    formfactor(i) = form_integral([1,1,1],250,theta(i));
end

figure(1)
clf; hold on;
plot(theta,formfactor)