%% b 5-40
R_halflife = 120;
P_halflife = 3600;
k_R = 0.01;
% b = 20;
b = 5:2:40;
eta = R_halflife/P_halflife;
gamma_p = log(2)/P_halflife;
gamma_r = gamma_p/eta;

p = k_R .* b./gamma_p;
fano = b./(1 + eta);
figure 
hold on
scatter(p, fano, 'ok', 'filled')
xlim([0 2000])
ylim([0 40])
xlabel('<p>')
ylabel('Fano factor ratio')
title('Controlling mean and variance of protein number')

%% k_R 0.0025-0.02
R_halflife = 120;
P_halflife = 3600;
%k_R = 0.01;
k_R = 0.0025:0.001:0.02;
b = 20;
eta = R_halflife/P_halflife;
gamma_p = log(2)/P_halflife;
gamma_r = gamma_p/eta;

p = k_R .* b./gamma_p;
fano = b./(1 + eta);
plot(p, fano, 'sk')

%% P_halflife 15min-2h
R_halflife = 120;
%P_halflife = 3600;
P_halflife = 900:200:7200;
k_R = 0.01;
b = 20;
eta = R_halflife./P_halflife;
gamma_p = log(2)./P_halflife;
gamma_r = gamma_p./eta;

p = k_R .* b./gamma_p;
fano = b./(1 + eta);
scatter(p, fano, 'dk', 'filled')

yline(1, '--');