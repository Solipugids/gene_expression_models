%% Monte carlo simulation
% transition: R --k_on--> A
% reverse transition: A--k_off--> R
% induced transcription: A --s_A--> M + A
% basal transcription: R --s_R--> M + R
% translation: M --s_P--> P + M
% mRNA decay: M --delta_M--> 0
% protein decay: P --delta_P--> 0
addpath('../Gillespie/')

%% Parameters
p.k_on = 10;
p.k_off = 10;
p.s_A = 50;
p.s_R = 5;
p.s_P = 0.2;
p.delta_M = 0.1;
p.delta_P = 0.05;
mol = 6.02214076e23; % mol-1
V = 200 * 1e-15; % m3
%% Initial state
tlast = [0, 20 * 60];%
x0    = [1, 0, 0, 0];% R, A, M, P

%% Specify reaction network
pfun = @propensities_2state;
stoich_matrix = [ -1 +1 0 0
                  +1 -1 0 0
                  0 0 +1 0
                  0 0 +1 0
                  0 0 0 +1
                  0 0 -1 0
                  0 0 0 -1];
 

%% Run simulation
[t,x] = directMethod(stoich_matrix, pfun, tlast, x0, p);
%[t_fig1a,x_fig1a] = firstReactionMethod(stoich_matrix, pfun, tlast, x0, p);

%% Plot time course
x_fig1a(:, 5) = (x_fig1a(:, 4)/mol/V) * 10e6;
figure();
hold on
stairs(t_fig1a,x_fig1a(:, 5)); 
set(gca,'XLim',tlast);
set(gca,'YLim',[0 0.15]);
xlabel('time (min)');
ylabel('Abundance(μM)');
title({'Low-amplitude fluctuations with high numbers of expressed mRNA', 'and protein molecules and a large cell volume.'})