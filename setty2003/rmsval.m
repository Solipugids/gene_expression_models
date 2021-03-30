function rms = rmsval(x, A, R, zdata)
    k_a = x(1);
    k_r = x(2);
    k_f = x(3);
    rms = sqrt(sum((zdata - (A.*(1-R).*k_a + R.*k_r + (1-A).*(1-R).*k_f)).^2));
%     v_1 = x(1);
%     v_2 = x(2);
%     v_3 = x(3);
%     v_4 = x(4);
%     v_5 = x(5);
%     rms = sqrt(sum(zdata - (v_1 .* (1 + v_2.*A + v_3.*R)./(1 + v_4.*A + v_5.*R))).^2);
end