%% parameter init
tlast = 20 * 60; % 20 hours totally. the unit of time is min.
dt = 1; % delta t
s_A = 50;
s_R = 5;
s_P = 0.2;
delta_M = 0.1;
delta_P = 0.05;
k_off = 10;
k_on = 10;
V = 200;


%% Deterministic
timepoints = 0:dt:tlast;
P_total = zeros(1, length(timepoints));
M_total = zeros(1, length(timepoints));

M_total(1) = 0;
P_total(1) = 0;
M = M_total(1);
P = P_total(1);
for i = 2: length(timepoints)
    dMdt = ((k_on/(k_off + k_on)) * (s_A/V)) + ((k_off/(k_off + k_on))*(s_R/V)) - (delta_M * M);
    dM = dMdt * dt;
    M_total(i) = M + dM;
    dPdt = s_P * M - delta_P * P;
    dP = dPdt * dt;
    P_total(i) = P + dP;
    M = M_total(i);
    P = P_total(i);
end
figure
plot(timepoints, P_total/10, 'LineWidth', 2)
ylim([0 1])