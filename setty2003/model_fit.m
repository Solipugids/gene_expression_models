camp = 0:0.1:10;
k_camp = 1.8;
n = 2;
A = zeros(length(camp), 1);
for c=1:length(camp)
    X = camp(c)/k_camp;
    A(c) = (X^n)/(1+(X^n));
end
h = gca;
plot(camp, A)
set(h, 'xscale', 'log')
xlabel("[cAMP]");
ylabel("Probability");

iptg = 0:1:100;
k_iptg = 1.2;
m = 4;
R = zeros(length(iptg),1);
for i=1:length(iptg)
    Y = iptg(i)/k_iptg;
    R(i) = 1/(1+(Y^m));
end
h = gca;
plot(iptg, R)
set(h, 'xscale', 'log')
xlabel("[IPTG]");
ylabel("Probability");

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
    pa(i, c) = f;
    s_pa = [s_pa; table(camp(c), iptg(i), pa(i,c))];
  end
end

camp = s_pa.Var1;
iptg = s_pa.Var2;
zdata = s_pa.Var3;
m = 4;
n = 2;
k_camp = 1.8;
k_iptg = 1.2;
X = camp./k_camp;
A = (X.^n)./(1+(X.^n));
Y = iptg./k_iptg;
R = 1./(1+(Y.^m));
type rmsval
fun = @(x)rmsval(x, A, R, zdata);
% x0 = [3.5, 70, 170, 17, 540];
x0 = rand(3, 1);
A = [];
b = [];
Aeq = [];
beq = [];
% lb = [0, 60, 140, 14, 440];
% ub = [5, 80, 200, 20, 640];
lb = [0 0 0];
ub = [100 100 100];
nonlcon = [];
% options = optimset('MaxFunEvals',100000, 'MaxIter', 10000);
options = optimoptions('fmincon','Algorithm','sqp','MaxIterations',1000);

bestx = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options)

%    k_a = 14.2313
%    k_r = 0
%    k_f = 8.9210
m = 4;
n = 2;
k_iptg = 1.2;
k_camp = 1.8;
iptg = 0:100;
camp = 0:0.1:10;
k_a = 14.2313;
k_r = 0;
k_f = 8.9210;
px = zeros(length(iptg), length(camp));
s_px = table();
for i=1:length(iptg)
  for c=1:length(camp)
    X = camp(c)/k_camp;
    A = (X^n)/(1+(X^n));
    Y = iptg(i)/k_iptg;
    R = 1/(1+(Y^m));
    f = (A.*(1-R).*k_a + R.*k_r + (1-A).*(1-R).*k_f);
    px(i, c) = log10(f);
    s_px = [s_px; table(c, i, px(i,c))];
  end
end


h = gca;
surf(px);
colormap jet;
shading interp;
set(h, 'xscale', 'log', 'yscale', 'log');
axis([0 200 0 200 0.5 1.5]);
xlabel("[cAMP]");
ylabel("[IPTG]");
zlabel("Promoter Activity");
xticks([2 11 100]);
xticklabels({'0.1', '1', '10'});
yticks([2 11 100]);
yticklabels({'1', '10', '100'});