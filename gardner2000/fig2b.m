figure
hold on

alpha_2 = 4;
gamma = 2;
u = 0:0.01:10;
v = alpha_2./(1 + u.^gamma);
plot(u, v, 'k', 'LineWidth', 2)

alpha_1 = 9;
beta = 2.1;
v = 0:0.01:10;
u = alpha_1./(1 + v.^beta);
plot(u, v, 'k', 'LineWidth', 2)

xlabel('u')
ylabel('v')
title('A bistable toggle network with balanced promoter strengths.')