%% Hejhej, Hanna.

theta = load("data/theta.tsv");
cross_section = load("data/cross_section.tsv");
meas_error = load("data/meas_error.tsv");

figure(1)
clf; hold on;
plot(theta,log10(cross_section))
