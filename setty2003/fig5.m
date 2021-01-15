m = 4;
n = 2;
k_iptg = 1.2;
k_camp = 1.8;
iptg = 0:100;
camp = 0:0.1:10;
V_1 = 3.5;
V_2 = 70;
V_3 = 170;
V_4 = 17;
V_5 = 540;
pa = zeros(length(iptg), length(camp));
s_pa = table();
for i=1:length(iptg)
  for c=1:length(camp)
    X = camp(c)/k_camp;
    A = (X^n)/(1+(X^n));
    Y = iptg(i)/k_iptg;
    R = 1/(1+(Y^m));
    f = V_1 * (1 + (V_2*A) + (V_3*R))/(1 + (V_4*A) + (V_5*R));
    pa(i, c) = log10(f);
    s_pa = [s_pa; table(c, i, pa(i,c))];
  end
end
s_pa.Properties.VariableNames = {'cAMP', 'IPTG', 'Activity'};
[smoothed, gof] = fit([s_pa.cAMP, s_pa.IPTG], s_pa.Activity, 'linearinterp');

h = gca;
surf(pa);
colormap jet;
shading interp;
set(h, 'xscale', 'log', 'yscale', 'log');
axis([0 200 0 200 0 1.5]);
xlabel("[cAMP]");
ylabel("[IPTG]");
zlabel("Promoter Activity");
xticks([2 11 100]);
xticklabels({'0.1', '1', '10'});
yticks([2 11 100]);
yticklabels({'1', '10', '100'});

h = gca;
pcolor(pa);
colormap jet;
shading interp;
set(h, 'xscale', 'log', 'yscale', 'log');
colorbar('AxisLocation','out', 'Ticks', [0.05 0.2 0.4 0.6 0.8 1], 'TickLabels', {'0', '0.2', '0.4', '0.6', '0.8', '1'});
xticks([2 11 100]);
xticklabels({'0.1', '1', '10'});
yticks([2 11 100]);
yticklabels({'1', '10', '100'});
xlabel("[cAMP]");
ylabel("[IPTG]");

%% smoothed
h = gca;
plot(smoothed);
colormap jet;
shading interp;
set(h, 'xscale', 'log', 'yscale', 'log');
axis([0 200 0 200 0 1.5]);
xlabel("[cAMP]");
ylabel("[IPTG]");
zlabel("Promoter Activity");
xticks([2 11 100]);
xticklabels({'0.1', '1', '10'});
yticks([2 11 100]);
yticklabels({'1', '10', '100'});

hold off
h = gca;
plot(smoothed, 'Style', 'Contour');
hold on
colormap jet;
shading interp;
set(h, 'xscale', 'log', 'yscale', 'log');
xticks([2 11 100]);
xticklabels({'0.1', '1', '10'});
yticks([2 11 100]);
yticklabels({'1', '10', '100'});
xlabel("[cAMP]");
ylabel("[IPTG]");
line([3.5 3.5], [9 101], "Color", "black", "LineWidth", 3);
line([3.5 20], [9 3.5], "Color", "black", "LineWidth", 3);
line([20 101], [3.5 3.5], "Color", "black", "LineWidth", 3);
line([1 3.5], [9 9], "Color", "black", "LineWidth", 3);
line([20 20], [3.5 1], "Color", "black", "LineWidth", 3);
hold off