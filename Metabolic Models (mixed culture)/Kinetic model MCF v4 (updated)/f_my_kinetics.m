function  r=f_my_kinetics(states,compounds,parameters, stoiMatrix, pH)

% This function calculates the reaction rates for each state
% *********************************************************************** %
% INPUTS:
% states 
% compounds 
% parameters
% stoichiometric matrix
% pH
% *********************************************************************** %
% OUTPUT:
% states reaction rate vector
% *********************************************************************** %

% Mateo Saavedra del Oso. CRETUS Institute. University of Santiago de Compostela.
% Spain. November 2020. Please contact msaavedra.deloso@usc.es if you intend 
% to use this code. 


statesNames=compounds.Abb;
kineticPar=parameters.kineticPar;
kineticNames=parameters.kineticNames;

H=10^(-pH);

%% States reaction rate vector 

r=0*states;


%% Compounds

Ssu = states(strcmp(statesNames,'Ssu'));        % Monosaccharides
Saa = states(strcmp(statesNames,'Saa'));        % Amino acids
Sfa = states(strcmp(statesNames,'Sfa'));        % LCFA
Sac = states(strcmp(statesNames,'Sac'));        % Acetate
Spro = states(strcmp(statesNames,'Spro'));      % Propionate
Sbu = states(strcmp(statesNames,'Sbu'));        % Butyrate
Sva = states(strcmp(statesNames,'Sva'));        % Valerate
Sh2 = states(strcmp(statesNames,'Sh2'));        % H2
Set = states(strcmp(statesNames,'Set'));        % Ethanol
Sch4 = states(strcmp(statesNames,'Sch4'));      % Methane
Stic = states(strcmp(statesNames,'Stic'));      % Total inorganic carbon
Sin = states(strcmp(statesNames,'Sin'));        % Inorganic nitrogen
Ssi = states(strcmp(statesNames,'Ssi'));        % Soluble inerts
Xc1 = states(strcmp(statesNames,'Xc1'));        % Composites 1
Xc2 = states(strcmp(statesNames,'Xc2'));        % Composites 2
Xc3 = states(strcmp(statesNames,'Xc3'));        % Composites from biomass decay
Xch = states(strcmp(statesNames,'Xch'));        % Carbohydrates
Xpr = states(strcmp(statesNames,'Xpr'));        % Proteins
Xtdf = states(strcmp(statesNames,'Xtdf'));      % Total dietary fiber: cellulose, hemicellulose, lignin
Xli = states(strcmp(statesNames,'Xli'));        % Lipids
Xi = states(strcmp(statesNames,'Xi'));          % Particulate inerts
Xaf = states(strcmp(statesNames,'Xaf'));        % Acidogenic fermentation biomass
Xfa = states(strcmp(statesNames,'Xfa'));        % LCFA degraders
Xc4 = states(strcmp(statesNames,'Xc4'));        % C4 degraders
Xpro = states(strcmp(statesNames,'Xpro'));      % Propionate degraders
Xac = states(strcmp(statesNames,'Xac'));        % Acetate degraders
Xh2 = states(strcmp(statesNames,'Xh2'));        % H2 degraders
Scat = states(strcmp(statesNames,'Scat'));      % Cations
San = states(strcmp(statesNames,'San'));        % Anions
Sgas_h2 = states(strcmp(statesNames,'Sgas_h2'));% H2 gas
Sgas_co2 = states(strcmp(statesNames,'Sgas_co2'));% H2 gas
Sgas_ch4 = states(strcmp(statesNames,'Sgas_ch4'));% H2 gas

%% Kinetic, Monod and inhibition parameters

kdis1 = kineticPar(strcmp(kineticNames,'kdis1'));         % Disintegration coeff. for composite 1
kdis2 = kineticPar(strcmp(kineticNames,'kdis2'));         % Disintegration coeff. for composite 2
kdis3 = kineticPar(strcmp(kineticNames,'kdis3'));         % Disintegration coeff. for composite 3
khyd_ch = kineticPar(strcmp(kineticNames,'khyd_ch'));     % Hydrolysis of carbohydrates coeff.
khyd_pr = kineticPar(strcmp(kineticNames,'khyd_pr'));     % Hydrolysis of proteins coeff.
khyd_li = kineticPar(strcmp(kineticNames,'khyd_li'));     % Hydrolysis of lipids coeff.
khyd_tdf = kineticPar(strcmp(kineticNames,'khyd_tdf'));   % Hydrolysis of TDF coeff.
km_su = kineticPar(strcmp(kineticNames,'km_su'));         % Maximum rate sugars acidogenic fermentation
km_aa = kineticPar(strcmp(kineticNames,'km_aa'));         % Maximum rate aminoacids acidogenic fermentation
km_fa = kineticPar(strcmp(kineticNames,'km_fa'));         % Maximum rate long fatty acids uptake
km_c4 = kineticPar(strcmp(kineticNames,'km_c4'));         % Maximum rate C4 uptade
km_pro = kineticPar(strcmp(kineticNames,'km_pro'));       % Maximum rate propionate uptake
km_ac = kineticPar(strcmp(kineticNames,'km_ac'));         % Maximum rate acetate uptake
km_h2 = kineticPar(strcmp(kineticNames,'km_h2'));         % Maximum rate hydrogen uptake
kdec_Xaf = kineticPar(strcmp(kineticNames,'kdec_Xaf'));   % Decay rate coeff. acidogenic fermentation biomass
kdec_Xfa = kineticPar(strcmp(kineticNames,'kdec_Xfa'));   % Decay rate coeff. long fatty degraders
kdec_Xc4 = kineticPar(strcmp(kineticNames,'kdec_Xc4'));   % Decay rate coeff. C4 degraders
kdec_Xpro = kineticPar(strcmp(kineticNames,'kdec_Xpro')); % Decay rate coeff. propionate degraders
kdec_Xac = kineticPar(strcmp(kineticNames,'kdec_Xac'));   % Decay rate coeff. acetate degraders
kdec_Xh2 = kineticPar(strcmp(kineticNames,'kdec_Xh2'));   % Decay rate coeff. hydrogen degraders
KS_su = kineticPar(strcmp(kineticNames,'KS_su'));         % Affinity constant sugars acidogenic fermentation
KS_aa = kineticPar(strcmp(kineticNames,'KS_aa'));         % Affinity constant aminoacids acidogenic fermentation
KS_fa = kineticPar(strcmp(kineticNames,'KS_fa'));         % Affinity constant long fatty acids uptake
KS_c4 = kineticPar(strcmp(kineticNames,'KS_c4'));         % Affinity constant C4 uptake
KS_pro = kineticPar(strcmp(kineticNames,'KS_pro'));       % Affinity constant propionate uptake
KS_ac = kineticPar(strcmp(kineticNames,'KS_ac'));         % Affinity constant acetate uptake
KS_h2 = kineticPar(strcmp(kineticNames,'KS_h2'));         % Affinity constant hydrogen uptake
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

%% Reactions rates

rho(1,1)=kdis1*Xc1;                                                              % Disintegration of composite 1 (rho1)
rho(2,1)=kdis2*Xc2;                                                              % Disintegration of composite 2 (rho2)
rho(3,1)=kdis3*Xc3;                                                              % Disintegration of composite from biomass decay (rho3)
rho(4,1)=khyd_ch*Xch;                                                            % Hydrolisis of carbohydrates (rho4)
rho(5,1)=khyd_pr*Xpr;                                                            % Hydrolisis of proteins (rho5)
rho(6,1)=khyd_tdf*Xtdf;                                                          % Hydrolisis of TDF (rho6)
rho(7,1)=khyd_li*Xli;                                                            % Hydrolisis of lipids (rho7)
rho(8,1)=km_aa*(Saa/(KS_aa+Saa))*Xaf*IpH_aa*Iin_lim;                             % Acidogenic fermentation of AA (rho8) 
rho(9,1)=km_su*(Ssu/(KS_su+Ssu))*Xaf*IpH_aa*Iin_lim;                             % Acidogenic fermentation of sugars (rho9) 
rho(10,1)=km_fa*(Sfa/(KS_fa+Sfa))*Xfa*IpH_aa*Iin_lim*Ih2_fa;                      % LCFA Uptake (rho10) 
rho(11,1)=km_c4*(Sva/(KS_c4+Sva))*Xc4*(Sva/(Sbu+Sva+1e-10))*IpH_aa*Iin_lim*Ih2_c4;% Valerate Uptake (rho11)
rho(12,1)=km_c4*(Sbu/(KS_c4+Sbu))*Xc4*(Sbu/(Sbu+Sva+1e-10))*IpH_aa*Iin_lim*Ih2_c4;% Butyrate Uptake (rho12)
rho(13,1)=km_pro*(Spro/(KS_pro+Spro))*Xpro*IpH_aa*Iin_lim*Ih2_pro;               % Propionaate Uptake (rho13) 
rho(14,1)=km_ac*(Sac/(KS_ac+Sac))*Xac*IpH_ac*Iin_lim*Inh3;                       % Acetate Uptake (rho14)
rho(15,1)=km_h2*(Sh2/(KS_h2+Sh2))*Xh2*IpH_h2*Iin_lim;                            % Hydrogen Uptake (rho15) 
rho(16,1)=kdec_Xaf*Xaf;                                                          % Decay of acidogenic fermentation biomass (rho16)
rho(17,1)=kdec_Xfa*Xfa;                                                          % Decay of LCFA (rho17)
rho(18,1)=kdec_Xc4*Xc4;                                                          % Decay of C4 degraders (rho18)
rho(19,1)=kdec_Xpro*Xpro;                                                        % Decay of propionate degraders (rho19)
rho(20,1)=kdec_Xac*Xac;                                                          % Decay of acetate degraders (rho20)
rho(21,1)=kdec_Xh2*Xh2;                                                          % Decay of hyrogen degraders (rho21)                          


%% Rate assigment to each state

r=stoiMatrix*rho;   % Update of states reaction rate vector

r(strcmp(statesNames,'Ssu'))=r(1);        % Monosaccharides
r(strcmp(statesNames,'Saa'))=r(2);        % Amino acids
r(strcmp(statesNames,'Sfa'))=r(3);        % LCFA
r(strcmp(statesNames,'Sac'))=r(4);        % Acetate
r(strcmp(statesNames,'Spro'))=r(5);       % Propionate
r(strcmp(statesNames,'Sbu'))=r(6);        % Butyrate
r(strcmp(statesNames,'Sva'))=r(7);        % Valerate
r(strcmp(statesNames,'Sh2'))=r(8);        % H2
r(strcmp(statesNames,'Set'))=r(9);        % Ethanol
r(strcmp(statesNames,'Sch4'))=r(10);      % Methane
r(strcmp(statesNames,'Stic'))=r(11);      % Total inorganic carbon
r(strcmp(statesNames,'Sin'))=r(12);       % Inorganic nitrogen
r(strcmp(statesNames,'Ssi'))=r(13);       % Soluble inerts
r(strcmp(statesNames,'Xc1'))=r(14);       % Composites 1
r(strcmp(statesNames,'Xc2'))=r(15);       % Composites 2
r(strcmp(statesNames,'Xc3'))=r(16);       % Composites from biomass decay
r(strcmp(statesNames,'Xch'))=r(17);       % Carbohydrates
r(strcmp(statesNames,'Xpr'))=r(18);       % Proteins
r(strcmp(statesNames,'Xtdf'))=r(19);      % Total dietary fiber: cellulose, hemicellulose, lignin
r(strcmp(statesNames,'Xli'))=r(20);       % Lipids
r(strcmp(statesNames,'Xi'))=r(21);        % Particulate inerts
r(strcmp(statesNames,'Xaf'))=r(22);       % Acidogenic fermentation biomass
r(strcmp(statesNames,'Xfa'))=r(23);       % LCFA degraders
r(strcmp(statesNames,'Xc4'))=r(24);       % C4 degraders
r(strcmp(statesNames,'Xpro'))=r(25);      % Propionate degraders
r(strcmp(statesNames,'Xac'))=r(26);       % Acetate degraders
r(strcmp(statesNames,'Xh2'))=r(27);       % H2 degraders
r(strcmp(statesNames,'Scat'))=r(28);      % Cations
r(strcmp(statesNames,'Xan'))=r(29);       % Anions
r(strcmp(statesNames,'Xgas_h2'))=r(30);   % Sgas_h2
r(strcmp(statesNames,'Xgas_co2'))=r(31);  % Sgas_co2
r(strcmp(statesNames,'Xgas_ch4'))=r(32);  % Sgas_ch4 

end
