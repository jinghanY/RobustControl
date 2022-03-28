function [A,B,G] = sysmatrix(Cf,Cr,V)
%%
lf = 1.1; 
lr = 1.58;
m = 1573;
Iz = 2873;
%%
a1 = -2 * (Cf + Cr)/(m * V);
a2 = 2 * (Cf + Cr)/m;
a3 = 2 * (-Cf * lf + Cr * lr)/(m * V);
a4 = -2 * (Cf * lf - Cr * lr)/(Iz * V);
a5 = 2 * (Cf * lf - Cr * lr)/(Iz);
a6 = -2 * (Cf * lf^2 + Cr * lr^2)/(Iz * V);
b1 = 2 * Cf /m;
b2 = 2 * Cf *lf / Iz;
g1 = -2 * (Cf * lf - Cr * lr)/(m * V) - V;
g2 = -2 * (Cf * lf^2 + Cr * lr^2)/(Iz * V);
%%
A = [0  1  0  0; 0 a1 a2 a3; 0  0  0  1; 0 a4 a5 a6 ];
B = [ 0; b1; 0; b2 ];
G = [ 0; g1; 0; g2 ];
end