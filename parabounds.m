%%
Omegaset=[Cf_lower/Cf_nominal Cf_upper/Cf_nominal];

%%
R_lower = 0.1;
SS = 1/(2*Cf_nominal * R_lower) * (2*Cf_upper*lf + Cr_upper*lr*(lr/lf -1)+ m*V^2/2);
Sigmaset=[-SS SS];
%%
XX_l = Cf_nominal/Cf_upper -1;
XX_u = Cf_nominal/Cf_lower -1;
XX = [XX_l XX_u];
%
TT1 = (km(1)*V * XX) /V ;
TT3 = -(km(3)*V * XX) /V ;
TT21 = (-((m+Iz)*(Cf_upper-Cf_nominal))/(2*m*Cf_lower) ...
       +((Iz*lr/lf-m)*(Cr_lower-Cr_nominal))/(2*m*Cf_upper) );

TT22 = (-((m+Iz)*(Cf_lower-Cf_nominal))/(2*m*Cf_upper) ...
       +((Iz*lr/lf-m)*(Cr_upper-Cr_nominal))/(2*m*Cf_lower));
TT2 = ([TT21 TT22] + km(2)*V * XX) /V ;

TT41 = -((m+Iz)*(Cf_upper-Cf_nominal))/(2*m*Cf_lower) ...
       +((Iz*lr/lf-m)*(Cr_lower-Cr_nominal))/(2*m*Cf_upper) ;
TT42 = -((m+Iz)*(Cf_lower-Cf_nominal))/(2*m*Cf_upper) ...
       +((Iz*lr/lf-m)*(Cr_upper-Cr_nominal))/(2*m*Cf_lower);
TT4 = ([TT41 TT42] + km(4)*V * XX) /V ;


Thetaset=[-1 1]*.2;   

