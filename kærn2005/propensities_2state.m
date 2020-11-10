function a = propensities_2state(x, p)
R = x(1);
A = x(2);
M = x(3);
P = x(4);

a = [p.k_on*R;
    p.k_off*A;
    p.s_A*A;
    p.s_R*R;
    p.s_P*M;
    p.delta_M*M;
    p.delta_P*P];
end