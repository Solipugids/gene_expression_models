t = 0:5000;
p_half = 3600;
b = 20;
k_R = 0.01;
eta = 0;
gamma_p = log(2)/p_half;

figure
hold on

fano = (( (1 - exp(-2 * gamma_p .* t))./(1 - exp(-gamma_p .* t)) ) .* b ./(1 + eta)) + 1;
plot(t, fano, 'k')

ylim([0 42])
xlabel('t(s)')
ylabel('fano')
title('Transient noise for a single gene.')

eta = 1/150;
fano = 1-exp(-gamma_p*t);
plot(t, fano)