function  Inhibition=f_my_inhibitions(Sin,Sh2,parameters, pH)

% This function calculates the reaction rates for each state
% *********************************************************************** %
% INPUTS:
% states 
% compounds name
% kinetic parameters
% stoichiometric matrix
% pH
% *********************************************************************** %
% OUTPUT:
% states reaction rate vector
% *********************************************************************** %

% Mateo Saavedra del Oso. CRETUS Institute. University of Santiago de Compostela.
% Spain. November 2020. Please contact msaavedra.deloso@usc.es if you intend 
% to use this code. 


kineticPar=parameters.kineticPar;
kineticNames=parameters.kineticNames;

H=10^(-pH);

%% Kinetic, Monod and inhibition parameters

KS_IN = kineticPar(strcmp(kineticNames,'KS_IN'));         % Affinity constant inorganic nitrogen
KI_NH3 = kineticPar(strcmp(kineticNames,'KI_NH3'));       % Ammonia inhibition
KI_h2_fa = kineticPar(strcmp(kineticNames,'KI_h2_fa'));   % Hydrogen inhibition for LFA uptake
KI_h2_c4 = kineticPar(strcmp(kineticNames,'KI_h2_c4'));   % Hydrogen inhibition for C4 uptake
KI_h2_pro = kineticPar(strcmp(kineticNames,'KI_h2_pro')); % Hydrogen inhibition for propionate uptake
pHUL_aa=kineticPar(strcmp(kineticNames,'pHUL_aa'));       % pH inhibition acidogenesis
pHLL_aa=kineticPar(strcmp(kineticNames,'pHLL_aa'));       % pH inhibition acidogenesis
pHUL_ac=kineticPar(strcmp(kineticNames,'pHUL_ac'));       % pH inhibition acetate uptake
pHLL_ac=kineticPar(strcmp(kineticNames,'pHLL_ac'));       % pH inhibition acetate uptake
pHUL_h2=kineticPar(strcmp(kineticNames,'pHUL_h2'));       % pH inhibition hydrogen uptake
pHLL_h2=kineticPar(strcmp(kineticNames,'pHLL_h2'));       % pH inhibition hydrogen uptake
Kac_in=kineticPar(strcmp(kineticNames,'Kac_in'));         % acidity constant acetate

%% Inhibition terms
Iin_lim = 1/(1+(KS_IN/(Sin+1e-10)));
Ih2_fa = 1/(1+(Sh2/KI_h2_fa)); 
Ih2_c4 = 1/(1+(Sh2/KI_h2_c4));
Ih2_pro = 1/(1+(Sh2/KI_h2_pro));
Snh3 = Sin*Kac_in/(H+Kac_in);
Inh3 = 1/(1+(Snh3/KI_NH3)); 

%pH Inhibition
naa = 3.0/(pHUL_aa-pHLL_aa);
nac = 3.0/(pHUL_ac-pHLL_ac);
nh2 = 3.0/(pHUL_h2-pHLL_h2);
KpH_aa=10^-(0.5*(pHUL_aa+pHLL_aa));
KpH_ac=10^-(0.5*(pHUL_ac+pHLL_ac));
KpH_h2=10^-(0.5*(pHUL_h2+pHLL_h2));
sH = 10^(-pH);



IpH_aa = KpH_aa^naa /(sH^naa + KpH_aa^naa);
IpH_ac = KpH_ac^nac /(sH^nac + KpH_ac^nac);
IpH_h2 = KpH_h2^nh2 /(sH^nh2 + KpH_h2^nh2);

Inhibition=[Iin_lim, Ih2_fa, Ih2_c4, Ih2_pro, Inh3, IpH_aa, IpH_ac, IpH_h2];