function [consumed_AA, consumed_glucose] = f_determine_my_consumed_AA_glucose(states, compounds, parameters,AA_profile, flagCase)

% This function determines the AA and glucose molar concentrations in the
% stationary state

% This function loads all parameters needed for the simulation of the VFA
% production from complex organic matter.
% *********************************************************************** %
% INPUTS:
% states
% compounds
% parameters
% AA_profile
% *********************************************************************** %
% OUTPUT:
% Aminoacids and glucose molar concentrations vector
% *********************************************************************** %

% Mateo Saavedra del Oso. CRETUS Institute. University of Santiago de Compostela.
% Spain. November 2020. Please contact msaavedra.deloso@usc.es if you intend 
% to use this code. 


%% Loading parameters and aminoacids composition

F=parameters.stoiPar;
FNames=parameters.stoiNames;
compoundsNames=compounds.Abb;
feeding=compounds.feed;
Fs_feed=compounds.Fs_feed;
Ffa_li=F(strcmp(FNames,'Ffa_li')); 

Fpr_xc1=F(strcmp(FNames,'Fpr_xc1'));    % protein fraction in composite 1
Fpr_xc2=F(strcmp(FNames,'Fpr_xc2'));    % protein fraction in composite 2
Fch_xc1=F(strcmp(FNames,'Fch_xc1'));    % carbohydrates fraction in composite 1
Fch_xc2=F(strcmp(FNames,'Fch_xc2'));    % carbohydrates fraction in composite 2
Ftdf_xc1=F(strcmp(FNames,'Ftdf_xc1'));  % TDF fraction in composite 1
Ftdf_xc2=F(strcmp(FNames,'Ftdf_xc2'));  % TDF fraction in composite 2
Fli_xc1=F(strcmp(FNames,'Fli_xc1'));  % TDF fraction in composite 1
Fli_xc2=F(strcmp(FNames,'Fli_xc2'));  % TDF fraction in composite 2

Xc1_feed=feeding(strcmp(compoundsNames,'Xc1'));  % Xc1 concentration in feed 
Xc2_feed=feeding(strcmp(compoundsNames,'Xc2'));  % Xc2 concentration in feed
Saa_feed=feeding(strcmp(compoundsNames,'Saa')); % Amino acids concentration in feed
Ssu_feed=feeding(strcmp(compoundsNames,'Ssu')); % Sugars concentration in feed

Xc1_ss=states(end,(strcmp(compoundsNames,'Xc1')));  % Xc1 concentration in reactor at stationary state
Xc2_ss=states(end,(strcmp(compoundsNames,'Xc2')));  % Xc2 concentration in reactor at stationary state
Xch_ss=states(end,(strcmp(compoundsNames,'Xch')));  % Xc1 concentration in reactor at stationary state
Xtdf_ss=states(end,(strcmp(compoundsNames,'Xtdf')));  % Xc1 concentration in reactor at stationary state
Xpr_ss=states(end,(strcmp(compoundsNames,'Xpr')));  % Xc2 concentration in reactor at stationary
Xli_ss=states(end,(strcmp(compoundsNames,'Xpr')));  % Xc2 concentration in reactor at stationary
Saa_ss=states(end,(strcmp(compoundsNames,'Saa')));  % Xc1 concentration in reactor at stationary state
Ssu_ss=states(end,(strcmp(compoundsNames,'Ssu')));  % Xc2 concentration in reactor at stationary

AA1=AA_profile(1,:);  % AA molar % composition for composite 1 
AA2=AA_profile(2,:);  % AA molar % composition for composite 2

MW_AA=AA_profile(4,:);   % AA molecular weight vector (g/mol)
COD_AA=AA_profile(5,:);  % AA chemical oxygen demand (g COD/L)


%% Consumed proteins and carbohydrates composition in stationary state

Ssu_eq_feed=Xc1_feed*(Fch_xc1+Ftdf_xc1+Fli_xc1*(1-Ffa_li))+Xc2_feed*(Fch_xc2+Ftdf_xc2+Fli_xc2*(1-Ffa_li))+Ssu_feed; %
Ssu_eq_ss=Xc1_ss*(Fch_xc1+Ftdf_xc1+Fli_xc1*(1-Ffa_li))+Xc2_ss*(Fch_xc2+Ftdf_xc2+Fli_xc2*(1-Ffa_li))+Xli_ss*(1-Ffa_li)+Xch_ss+Xtdf_ss+Ssu_ss; % 
Ssu_consumed=Ssu_eq_feed-Ssu_eq_ss;

if Ssu_consumed<0
    Ssu_consumed=0;
end

Xpr_xc1_cons=(Xc1_feed-Xc1_ss)*Fpr_xc1;
Xpr_xc2_cons=(Xc2_feed-Xc2_ss)*Fpr_xc2;

if flagCase.simulation==2
    Saa_xc1_cons=sum(feeding)*(1-flagCase.substrates_ratioj)*Fpr_xc1*Fs_feed(flagCase.substrate1);
    Saa_xc2_cons=sum(feeding)*flagCase.substrates_ratioj*Fpr_xc2*Fs_feed(flagCase.substrate2);
else
    Saa_xc1_cons=sum(feeding)*(1-flagCase.substrates_ratio)*Fpr_xc1*Fs_feed(flagCase.substrate1);
    Saa_xc2_cons=sum(feeding)*flagCase.substrates_ratio*Fpr_xc2*Fs_feed(flagCase.substrate2);
end

Fsaa_xc1=(Xpr_xc1_cons+Saa_xc1_cons)/(Xpr_xc1_cons+Xpr_xc2_cons+Saa_xc1_cons+Saa_xc2_cons);

Saa_consumed=Xpr_xc1_cons+Xpr_xc2_cons+Saa_xc1_cons+Saa_xc2_cons-Xpr_ss-Saa_ss;

if Saa_consumed<0
    Saa_consumed=0;
end
consumed_AA=Saa_consumed*(Fsaa_xc1*AA1.*COD_AA/sum(AA1.*COD_AA)+(1-Fsaa_xc1)*AA2.*COD_AA/sum(AA2.*COD_AA))./COD_AA;

% Consumed AA vector 
consumed_AA=consumed_AA';
% Consumed AA molar concentration vector 

consumed_glucose=Ssu_consumed/192; % consumed glucose molar concentration 

end

