function [V_opt] = V_max(Cf_lower, Cf_upper, Cr_lower, Cr_upper)
    Cf = 31234;
    Cr = 31234;

    Cf_nominal = (Cf_lower + Cf_upper)/2;
    Cr_nominal = (Cr_lower + Cr_upper)/2;
    %
    R = 20; % meter
    V = 30; % 30 m/s max
    m = 1573;
    Iz = 2873;
    lf = 1.1; 
    lr = 1.58;
    %
    km = [0.7223; 2.5855; -0.6669; 0.1873];
    [A0,bm,g0] = sysmatrix(Cf_nominal,Cr_nominal,V);
    Am = A0 - bm*km';
    % Unknown parameter bounds
    %%%%%%%%%%%%%
    parabounds;
    %%%%%%%%%%%%%
    P  =lyap(Am',eye(4));
    Gamma = 10000;
    k = 20;
    %
    vOpt
end 


