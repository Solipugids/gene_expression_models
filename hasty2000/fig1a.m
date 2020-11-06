x = 0:0.01:0.8;
alpha = 50;
gammaMax = 20;
y = (alpha .* (x.^2))./(1 + (2 .* (x.^2)) + (5 .* (x.^4)));
figure
hold on
plot(x, y, 'k', 'LineWidth', 2)
for gamma = 10:2:gammaMax
    z = (gamma .* x) -1;
    plot(x, z, 'b--')
end
xlabel('x')
ylabel('f(x)')
title('Graphical depiction of the fixed points of Eq. by x = 0')