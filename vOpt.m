Opara.TT1 = TT1;
Opara.TT2 = TT2;
Opara.TT3 = TT3;
Opara.TT4 = TT4;
Opara.Cf = Cf;
Opara.Cr = Cr;
Opara.Cf_nominal = Cf_nominal;
Opara.Cr_nominal = Cr_nominal;
Opara.Cf_lower = Cf_lower;
Opara.Cf_upper = Cf_upper;
Opara.Cr_lower = Cr_lower;
Opara.Cr_upper = Cr_upper;
Opara.R_lower = R_lower;
Opara.P = P;
Opara.Iz = Iz;
Opara.m = m;
Opara.lf = lf;
Opara.lr = lr;
Opara.Gamma = Gamma;
%bounds
Opara.b1= 1.1;
Opara.b2= 1.1;
Opara.b3= 1.1;
%%
clc
A = []; b = []; Aeq = []; beq = [];
lb = [0,  TT1(1), TT2(1), TT3(1), TT4(1), Omegaset(1)];
ub = [30, TT1(2), TT2(2), TT3(2), TT4(2), Omegaset(2)];
x0 = [10, 0.01, 0.01, 0.01, 0.01, 1];
x00 = [30, 0.01, 0.01, 0.01, 0.01, 1];
%%  HL1
fun1 = @(x) HGH1_norm(x,Opara);
x1 = fmincon(fun1,x0,A,b,Aeq,beq,lb,ub);
HL1 = HGH1_norm(x1,Opara); % 2
%% CL1
fun2 = @(x) C_norm(x,Opara);
x2 = fmincon(fun2,x0,A,b,Aeq,beq,lb,ub);
CL1 = C_norm(x2,Opara); % 0
%% H_L1
fun3 = @(x) H_L1_norm(x,Opara);
x3 = fmincon(fun3,x0,A,b,Aeq,beq,lb,ub)
H_L1 = H_L1_norm(x3,Opara); % 48.9684
%%
Opt.HL1 = HL1;
Opt.CL1 = CL1;
Opt.H_L1 = H_L1;
%%
clc
fun4 = @(V) -V;
lbv = [1];
ubv = [30];
V_opt = fmincon(fun4,15,A,b,Aeq,beq,lbv,ubv,@(V)V_con(V,Opara,Opt));
%%
% HGH1_norm(x0,Opara)
%%
function [HL1_norm] = HGH1_norm(x,Opara)
    Cf_nominal = Opara.Cf_nominal;
    Cr_nominal = Opara.Cr_nominal;
    V = x(1);
    theta = [x(2) x(3) x(4) x(5)];
    omega = x(6);
    %
    [A0,bm,g0] = sysmatrix(Cf_nominal,Cr_nominal,V);
     km = [0.7223; 2.5855; -0.6669; 0.1873];
     Am = A0 - bm*km';
     k = 20;
     Ag = [Am + bm*theta bm*omega; -k*theta -k*omega];

     s = tf('s');
     Hg = inv(s*eye(5) - Ag) * [bm;0];
     %
     Cs = (omega*k/s)/(1+omega*k/s);
     Hs = inv(s*eye(4) - Am) * bm;
     H1 = Cs/(ones(1,4)*Hs) * ones(1,4);
     %
     HL1_norm = norm(Hg*H1,Inf,1e-10);
end
%%
function [CL1_norm] = C_norm(x,Opara)
    Cf_nominal = Opara.Cf_nominal;
    Cr_nominal = Opara.Cr_nominal;
    V = x(1);
    theta = [x(2) x(3) x(4) x(5)];
    omega = x(6);
    %
    [A0,bm,g0] = sysmatrix(Cf_nominal,Cr_nominal,V);
     km = [0.7223; 2.5855; -0.6669; 0.1873];
     Am = A0 - bm*km';
     k = 20;
     Ag = [Am + bm*theta bm*omega; -k*theta -k*omega];

     s = tf('s');
     Cs = (omega*k/s)/(1+omega*k/s);

     CL1_norm = norm(Cs / omega * theta,Inf,1e-10);

end
%%
function [H_L1_norm] = H_L1_norm(x,Opara)
    Cf_nominal = Opara.Cf_nominal;
    Cr_nominal = Opara.Cr_nominal;
    V = x(1);
    omega = x(2);
    %
     k = 20;
     [A0,bm,g0] = sysmatrix(Cf_nominal,Cr_nominal,V);
     km = [0.7223; 2.5855; -0.6669; 0.1873];
     Am = A0 - bm*km';

     s = tf('s');
     %
     Hs = inv(s*eye(4) - Am) * bm;
     Cs = (omega*k/s)/(1+omega*k/s);
     H1 = Cs/(ones(1,4)*Hs) * ones(1,4);
     %
     H_L1_norm = norm(H1/omega,Inf,1e-10);
end

%%
function [c,ceq] = V_con(V,Opara,Opt)
    %%
    TT1 = Opara.TT1;TT2 = Opara.TT2;TT3 = Opara.TT3;TT4 = Opara.TT4;
    Cf_upper = Opara.Cf_upper;Cf_lower = Opara.Cf_lower;
    Cr_upper = Opara.Cr_upper;Cr_lower = Opara.Cr_lower;
    Cf_nominal = Opara.Cf_nominal;Cr_nominal = Opara.Cr_nominal;
    R_lower = Opara.R_lower;
    P = Opara.P;
    m = Opara.m;
    lf = Opara.lf;
    lr = Opara.lr;
    Gamma = Opara.Gamma;
    b1 = Opara.b1;
    b2 = Opara.b2;
    b3 = Opara.b3;
    HL1 = Opt.HL1;
    CL1 = Opt.CL1;
    H_L1 = Opt.H_L1;
    %%
    Vmin = 1;
    Vmax = 30;
    %%
    [A0,bm,g0] = sysmatrix(Cf_nominal,Cr_nominal,Vmin);
     km = [0.7223; 2.5855; -0.6669; 0.1873];
     Am_min = A0 - bm*km';
    [A0,bm,g0] = sysmatrix(Cf_nominal,Cr_nominal,Vmax);
     Am_max = A0 - bm*km';

    Qmin = -(Am_min * P + P' * Am_min);
    Qmax = -(Am_max * P + P' * Am_max);
    Rd = 40;
    minQmin = min(abs(eig(Qmin)));
    minQmax = min(abs(eig(Qmax)));
    %%
    c = 4/V^2 * ((TT1(2))^2 + (TT2(2))^2 + (TT3(2))^2 + (TT4(2))^2) ...
        + ((Cf_upper - Cf_lower)/Cf_nominal)^2 ...
        + 1/((Cf_nominal)^2 - R_lower) ...
        *( min(abs(eig(P))) * Rd / (min(minQmin,minQmax)) + 1/R_lower)...
        * (2*Cf_upper*lf + Cr_upper*lr*(lr/lf-1)+m*V^2)^2 ...
        - max(abs(eig(P))) ...
            * Gamma * min([b1^2, b2^2/(HL1),b3^2/(CL1*HL1+H_L1)]);
    %%
    ceq =[];
end
