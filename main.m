%% Constants
global hbar q_e m_e c alpha eV millibarn giga mega fermi
% Unit conversion
eV        = 1.602e-19;
millibarn = 1e-31;
mega      = 1e6;
giga      = 1e9;
fermi     = 1e-15;

% Fundamental constants 
% All taken from Nuclear and Particle physics 2nd edition
c     = 2.998e8;
hbar  = 6.582e-25*giga*eV;
q_e   = 1.602e-19;
m_e   = 0.511*mega*eV/c^2;
alpha = 1/(137.04);

% General relativistic functions
gamma   = @(E) E/(m_e*c^2);
p_e     = @(E) sqrt((E/c).^2 - m_e^2*c^2); 
v       = @(E) p_e(E)./(m_e*gamma(E));
beta    = @(E) v(E)/c;

% Measurement data
theta         = load("data/theta.tsv");
theta_rad     = theta*pi/180;
cross_section = load("data/cross_section.tsv");
meas_error    = load("data/meas_error.tsv");

% Experimental values 
Z_Ca40   = 20;          % Number charge of Calcium nucleus
Z_e      = 1;           % Number charge of electron
E_250MeV = 250*mega*eV; % Energy in electron beam

% Vector used for passing data to functions
data       = [theta_rad, cross_section, meas_error];
exp_values = [Z_Ca40, Z_e, E_250MeV];

% Upper limit us for integration
% It is possible to set a fixed upper limit for all integrals in this
% assignment since the all have a exponentially decaying term
upper_limit = 100e-15;

% Making pretty figures
In = 'Interpreter';
La = 'Latex';

%% Proton density and fitting vector
% Create initial guesses for the fitting!
% rho_0 highest proton density
% a middle of fall off
% b sharpness of the falls off
rho_0   = q_e*0.07*1e45; 
a_0     = 3e-15;  
b_0     = 0.5e-15;

% Scaling the optimazation parameters. The paramters for rho_0 and a,b are
% approx 10^(40) apart. We rescale them so lsqcurvefit has to optimize
% values that are of equal size.
X_real  = [rho_0, a_0, b_0];
X_scale = [1/rho_0, 1/a_0, 1/b_0];
X_0     = X_real.*X_scale;

% Proton density that we want to optimize
% X(1) = rho_0, X(2) = a, X(3) = b
rho_ch  = @(X, r) (X(1)/X_scale(1))./(1 + exp((r-(X(2)/X_scale(2)))./(X(3)/X_scale(3))));

% Plotting the proton desnsity for the initial guess, X_0
r = linspace(0,8,1000)*fermi;

figure(1)
clf; hold on;
e_fm3 = q_e*1e45; 
plot(r/fermi, rho_ch(X_0, r)/e_fm3,'blue','Linewidth',2)

title('Inital guess for $\rho_{ch}$',In,La,'Fontsize',22)
xlabel('$r$ [fm]',In,La,'Fontsize',18)
ylabel('$\rho$ [e$\textrm{fm}^{-3}$]',In,La,'Fontsize',18)
legend('$\rho_{ch}$',In,La,'Fontsize',20)

%% Theoretical cross-section
% Difference of p_e before and after collision with Ca40
% Assumes the heavy core to be stationary and not recoil from collision 
q_diff = @(E, theta) sqrt(4*p_e(E).^2*sin(theta/2).^2);

% The form factor
form_const      = @(E, theta) 4*pi*hbar./(Z_Ca40*q_e*q_diff(E, theta));
form_integrand  = @(X, r, E, theta) r.*rho_ch(X, r).*sin(q_diff(E, theta).*r/hbar);
form_integral   = @(X, E, theta) integral(@(r) form_integrand(X, r, E, theta), 0, upper_limit,'ArrayValued',true);
formfactor      = @(X, E, theta) form_integral(X, E, theta).*form_const(E, theta);

% The Rutherford and Mott cross-sections
cross_ruth      = @(E, theta) (Z_e*Z_Ca40*alpha*hbar*c)^2./(4*beta(E)^4*E^2*sin(theta/2).^4);
cross_mott      = @(E, theta) cross_ruth(E, theta).*(1 - beta(E)^2*sin(theta/2).^2);

% Together they make the theoretical cross-section to be fitted to the
% measurement
cross_theo      = @(X, E, theta) cross_mott(E, theta).*abs(formfactor(X, E, theta)).^2;

%% Optimze the theoretical cross-section to the experiment
% x-data for the fit. The extra added zero is since the charge
% conservation in independet of theta
xdata = [theta_rad;0];

% Weight for the conservation of charge
eps_charge = 1e-6;

% y-data for the fit. All of the measurements plus the weighted charge.
ydata = [cross_section./meas_error; Z_Ca40/eps_charge];

% Options for the lsqcurvefit. The options are to use the
% levenberg-marquardt since it seemed to converge faster for our problem and 
% display information when interating. 
% The rest of the settings is to make the optimizer work harder.
options_lsq = optimoptions('lsqcurvefit','Algorithm',    'levenberg-marquardt',...
                                         'Display',      'iter',...
                                         'MaxFunEvals',  1e8,...
                                         'MaxIter',      1e6,...
                                         'TolFun',       1e-10,...
                                         'TolX',         1e-20);

% Optimizes optimize_lsq for X_0.
X_star = lsqcurvefit(@(X,xdata) optimize_lsq(X, xdata, rho_ch, cross_theo, data, exp_values,...
                                eps_charge,...
                                upper_limit,...
                                millibarn, q_e),...
                                X_0,...
                                xdata,...
                                ydata,...
                                [],[],...
                                options_lsq);
                           
X_star_real = X_star./X_scale;  % Scale the optimzied parameters to find the "real" ones

% Plot and check the theoretic cross-section with X_star
% Calculate the found rms charge
integrand       = @(X, r) r.^4.*rho_ch(X, r);
radius_squared  = (4*pi/(exp_values(1)*q_e))*integral(@(r) integrand(X_star, r), 0, upper_limit);
radius_X_star   = sqrt(radius_squared);

% Calculate the total charge with X_star 
integrand       = @(X, r) r.^2.*rho_ch(X, r);
total_charge    = @(X) (4*pi/q_e)*integral(@(r) integrand(X, r), 0, upper_limit);
charge_X_star   = total_charge(X_star);

% Plot the cross-section with the optimized charge distriutions vs the 
% cross-section from the experiment. Figure 2.1 in the assignment
theta_cont = linspace(theta(1)-10,theta(end)+20,1000);
theta_rad_cont = theta_cont*pi/180;
figure(2)
clf; hold on;
plot(theta, log10(cross_section),'redo','MarkerFaceColor','r','MarkerSize',6)
plot(theta_cont, log10(cross_theo(X_star, E_250MeV,theta_rad_cont)/millibarn),'black','linewidth',2)
axis([20,140,-8,0])

title('$\sigma_e$ and $\sigma_{X^*}$',In,La,'Fontsize',22)
xlabel('$\theta$ [deg]',In,La,'Fontsize',18)
ylabel('$d\sigma$/$d\Omega$ [mb/sr]',In,La,'Fontsize',18)
legend('$\sigma_e$','$\sigma_{X^*}$',In,La,'Fontsize',20)

% Plot the optimized charge distribution. Figure 2.2 in the assignment
e_fm3 = q_e*1e45; 
figure(3)
clf; hold on;
plot(r/fermi, rho_ch(X_star, r)/e_fm3,'blue','linewidth',2)

title('Optimized charge distribution $\rho_{ch}(r;X_*)$',In,La,'Fontsize',22)
xlabel('$r$ [fm]',In,La,'Fontsize',18)
ylabel('$\rho$ [e$\textrm{fm}^{-3}$]',In,La,'Fontsize',18)
legend('$\rho$',In,La,'Fontsize',20)

% Pretty output
fprintf('RMS charge radius\n\n')
fprintf('------ r:\t %1.4f fm\n\n', radius_X_star/fermi)

fprintf('Charge check\n\n')
fprintf('------ Z:\t %1.12f\n\n', charge_X_star)

fprintf('Optimized parameters\n\n')
fprintf('------ rho_ch:\t %1.4d\n', X_star_real(1))
fprintf('------ a:\t %1.4d\n', X_star_real(2))
fprintf('------ b:\t %1.4d\n', X_star_real(3))

