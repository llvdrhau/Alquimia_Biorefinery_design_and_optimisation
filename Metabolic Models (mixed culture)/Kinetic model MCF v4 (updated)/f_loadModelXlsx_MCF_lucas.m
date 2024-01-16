function [compounds, parameters, stoiMatrix, reaction, AA_profile,...
    newParticulateSoluble, newCarbohydrateProteinScharbohydrates, newCompositionAA] = ...
    f_loadModelXlsx_MCF_lucas(flagCase,i, feed_concentration)

% This function loads all parameters needed for the simulation of the VFA
% production from complex organic matter.
% *********************************************************************** %
% INPUTS:
% myMCFModelStructure.xlsx
% *********************************************************************** %
% OUTPUT:
% compounds name and index, and initial concentrations
% feeding 
% kinetic paremeters
% stoichiometric parameters
% stoichiometric matrix
% *********************************************************************** %

% Mateo Saavedra del Oso. CRETUS Institute. University of Santiago de Compostela.
% Spain. November 2020. Please contact msaavedra.deloso@usc.es if you intend 
% to use this code. 

flagCase=flagCase;
warning off all

route = char('myMCFModelStructure.xlsx'); 


%% load states Abb, index, initial conditions, feeding concentration, MW and flags

[data, names]=xlsread(route, strcat('Compounds'));  % read states index, names and initial conditions

index=data(:,1);                                  % compounds index
Abb=names(2:end,1);                               % states names
initial=data(:,2);                                % states initial conditions
feed=data(:,3);                                   % feeding states concentration
MW=data(:,4);                                     % molecular weight
flagLiquid=logical(data(:,5));                    % flag liquid
flagGas=logical(data(:,6));                       % flag gas


%% load parameters

[par,parNames]= xlsread(route, strcat('Parameters')); 
% read kinetic parameters from a xls file

% Reactor parameters and constants
aux = strcmp('Reactor',parNames(1,:));
reactorPar=par((isfinite(par(:,aux))),aux); 
reactorNames=parNames(2:end,aux);                   
reactorNames(strcmp('',reactorNames)) = [];           

% Kinetic parameters
aux = strcmp('Kinetics',parNames(1,:));
kineticPar=par((isfinite(par(:,aux))),aux); 
kineticNames=parNames(2:end,aux);                 
kineticNames(strcmp('',kineticNames)) = [];           

% Stoichiometric parameters
aux = strcmp('Stoichiometric',parNames(1,:));
stoiPar=par((isfinite(par(:,aux))),aux); 
stoiNames=parNames(2:end,aux);                 
stoiNames(strcmp('',stoiNames)) = [];   


%% load substrates library

[parSubstrate, substrateName]=xlsread(route, strcat('Feed'));
% read substrates library

substrateNames=substrateName(1,3:end);
parametersNames=substrateName(3:end,1);
substrateIndex=parSubstrate(1,2:end);
parameterIndex=parSubstrate(2:end,1);
parSubstrate=parSubstrate(2:end,:);
pos_AAini=parameterIndex(strcmp(parametersNames,'AAini'));
pos_AAend=parameterIndex(strcmp(parametersNames,'AAend'));
pos_Fxc_feed=parameterIndex(strcmp(parametersNames,'Fxc'));
pos_Fx_feed=parameterIndex(strcmp(parametersNames,'Fx'));
pos_Fs_feed=parameterIndex(strcmp(parametersNames,'Fs'));
Fxc_feed=parSubstrate(pos_Fxc_feed,2:end);
Fx_feed=parSubstrate(pos_Fx_feed,2:end);
Fs_feed=parSubstrate(pos_Fs_feed,2:end);
composition_matrix=parSubstrate((pos_Fs_feed+1):(pos_AAini-1),2:end);
AA_pool=parSubstrate(pos_AAini:pos_AAend,2:end);
kinetics=parSubstrate(end,2:end);

stoiPar(1:8)=composition_matrix(:,flagCase.substrate1); %assign composition of s1 to stoiPar
stoiPar(9:16)=composition_matrix(:,flagCase.substrate2); %assign composition of s2 to stoiPar



%% update the feeding

if length(flagCase.substrates_ratio)==1
    Xc1=feed_concentration*(1-flagCase.substrates_ratio)*Fxc_feed(flagCase.substrate1);
    Xc2=feed_concentration*flagCase.substrates_ratio*Fxc_feed(flagCase.substrate2);
    X_macromolecules=feed_concentration*((1-flagCase.substrates_ratio)*Fx_feed(flagCase.substrate1).*composition_matrix(1:4,flagCase.substrate1)+flagCase.substrates_ratio*Fx_feed(flagCase.substrate2).*composition_matrix(1:4,flagCase.substrate2));
    aux=zeros(3,length(substrateNames));
    aux(1,:)=composition_matrix(1,:)+composition_matrix(3,:);
    aux(2,:)=composition_matrix(2,:);
    aux(3,:)=composition_matrix(4,:);
    S_monomers=feed_concentration*((1-flagCase.substrates_ratio)*Fs_feed(flagCase.substrate1).*aux(:,flagCase.substrate1)+flagCase.substrates_ratio*Fs_feed(flagCase.substrate2).*aux(:,flagCase.substrate2));
    feed(strcmp(Abb,'Xc1'))=Xc1;
    feed(strcmp(Abb,'Xc2'))=Xc2;
    feed(strcmp(Abb,'Ssu'))=S_monomers(1);
    feed(strcmp(Abb,'Saa'))=S_monomers(2);
    feed(strcmp(Abb,'Sfa'))=S_monomers(3);
    feed(strcmp(Abb,'Xch'))=X_macromolecules(1);
    feed(strcmp(Abb,'Xpr'))=X_macromolecules(2);
    feed(strcmp(Abb,'Xtdf'))=X_macromolecules(3);
    feed(strcmp(Abb,'Xli'))=X_macromolecules(4);
    
    Cxc1=stoiPar(strcmp(stoiNames,'Cxc1'));
    Cxc2=stoiPar(strcmp(stoiNames,'Cxc2'));
    Nxc1=stoiPar(strcmp(stoiNames,'Nxc1'));
    Nxc2=stoiPar(strcmp(stoiNames,'Nxc2'));
    
    Stic=Xc1*Cxc1+Xc2*Cxc2;
    Sin=Xc1*Nxc1+Xc2*Nxc2;
    
    feed(strcmp(Abb,'Stic'))=Stic;
    feed(strcmp(Abb,'Sin'))=Sin;
    
    
    
    
    
else
    Xc1=feed_concentration*(1-flagCase.substrates_ratio(i))*Fxc_feed(flagCase.substrate1);
    Xc2=feed_concentration*flagCase.substrates_ratio(i)*Fxc_feed(flagCase.substrate2);
    X_macromolecules=feed_concentration*((1-flagCase.substrates_ratio(i))*Fx_feed(flagCase.substrate1).*composition_matrix(1:4,flagCase.substrate1)+flagCase.substrates_ratio(i)*Fx_feed(flagCase.substrate2).*composition_matrix(1:4,flagCase.substrate2));
    aux=zeros(3,length(substrateNames));
    aux(1,:)=composition_matrix(1,:)+composition_matrix(3,:);
    aux(2,:)=composition_matrix(2,:);
    aux(3,:)=composition_matrix(4,:);
    S_monomers=feed_concentration*((1-flagCase.substrates_ratio(i))*Fs_feed(flagCase.substrate1).*aux(:,flagCase.substrate1)+flagCase.substrates_ratio(i)*Fs_feed(flagCase.substrate2).*aux(:,flagCase.substrate2));
    feed(strcmp(Abb,'Xc1'))=Xc1;
    feed(strcmp(Abb,'Xc2'))=Xc2;
    feed(strcmp(Abb,'Ssu'))=S_monomers(1);
    feed(strcmp(Abb,'Saa'))=S_monomers(2);
    feed(strcmp(Abb,'Sfa'))=S_monomers(3);
    feed(strcmp(Abb,'Xch'))=X_macromolecules(1);
    feed(strcmp(Abb,'Xpr'))=X_macromolecules(2);
    feed(strcmp(Abb,'Xtdf'))=X_macromolecules(3);
    feed(strcmp(Abb,'Xli'))=X_macromolecules(4);
    
    Cxc1=stoiPar(strcmp(stoiNames,'Cxc1'));
    Cxc2=stoiPar(strcmp(stoiNames,'Cxc2'));
    Nxc1=stoiPar(strcmp(stoiNames,'Nxc1'));
    Nxc2=stoiPar(strcmp(stoiNames,'Nxc2'));
    
    Stic=Xc1*Cxc1+Xc2*Cxc2;
    Sin=Xc1*Nxc1+Xc2*Nxc2;
    
    feed(strcmp(Abb,'Stic'))=Stic;
    feed(strcmp(Abb,'Sin'))=Sin;
end


%% load stoichiometric table and build stoichiometric matrix

[~, ~, stoiTable]=xlsread(route, strcat('ReactionMatrixFull'));
% read stoichiometry reaction matrix

reaction.Names=stoiTable(2,3:end-1);
reaction.Index=stoiTable(1,3:end-1);

stoiMatrix = f_my_stoichiometric_matrix(stoiTable, stoiPar, stoiNames);
% this function builds the stoichiometric matrix


%% read AA profile

route = char('myMCFModelStructure.xlsx');
AA_profile=xlsread(route, strcat('AA profile'));  % read AA composition for each composite
AA_profile(1,:)=AA_pool(:,flagCase.substrate1)';
AA_profile(2,:)=AA_pool(:,flagCase.substrate2)';

% AA in that order:
% Arg Ala Asp Lys Glut Ser Thr Cys Gly Pro Val IsoL Leu Meth GluM AspG Hist

%% function output

compounds.index=index;
compounds.Abb=Abb;
compounds.initial=initial;
compounds.feed=feed;
compounds.MW=MW;
compounds.flagLiquid=flagLiquid;
compounds.flagGas=flagGas;
compounds.Fs_feed=Fs_feed;

parameters.reactorPar=reactorPar;
parameters.reactorNames=reactorNames;
parameters.kineticPar=kineticPar;
parameters.kineticNames=kineticNames;
parameters.stoiPar=stoiPar;
parameters.stoiNames=stoiNames;


%% calculate new fractions for
% Added by Lucas Van der Hauwaert
% lucas.vanderhauwaert@usc.es 
% determines combined different fractions of the two streams 
% fraction of 1) charbohydrates/shortChainCharbohydrates/protein/ 
%             2) particulates/soluble
%             3) protein composition 

% particular - soluble 
composition_particulate_soluble_stream1 = [parSubstrate(strcmp(parametersNames,'Fxc'),flagCase.substrate1+1), parSubstrate(strcmp(parametersNames,'Fs'),flagCase.substrate1+1)];
composition_particulate_soluble_stream2 = [parSubstrate(strcmp(parametersNames,'Fxc'),flagCase.substrate2+1), parSubstrate(strcmp(parametersNames,'Fs'),flagCase.substrate2+1)];

% carbohydrates - proteins - shortCarbs
composition_carb_protein_sCarb_stream1 = [parSubstrate(strcmp(parametersNames,'Fch'),flagCase.substrate1+1),... % carbohydrates
    parSubstrate(strcmp(parametersNames,'Fpr'),flagCase.substrate1+1), ... % protein
    parSubstrate(strcmp(parametersNames,'Ftdf'),flagCase.substrate1+1)]; % quickly degraded carbohydrates

composition_carb_protein_sCarb_stream2 = [parSubstrate(strcmp(parametersNames,'Fch'),flagCase.substrate2+1), ... % carbohydrates
    parSubstrate(strcmp(parametersNames,'Fpr'),flagCase.substrate2+1), ... % protein
    parSubstrate(strcmp(parametersNames,'Ftdf'),flagCase.substrate2+1)]; % quickly degraded carbohydrates
% protein compositions 
composition_1_AA = AA_profile(1,:);
composition_2_AA = AA_profile(2,:);

% recalucualte new compostions 
splitRatio = flagCase.substrates_ratio;
FComposition1 = splitRatio*feed_concentration;
FComposition2 = (1-splitRatio)*feed_concentration;

newParticulateSoluble = (FComposition1.*composition_particulate_soluble_stream1 + FComposition2.*composition_particulate_soluble_stream2)./ feed_concentration;
newCarbohydrateProteinScharbohydrates = (FComposition1.*composition_carb_protein_sCarb_stream1 + FComposition2.*composition_carb_protein_sCarb_stream2)./ feed_concentration;
newCompositionAA = (FComposition1.*composition_1_AA + FComposition2.*composition_2_AA)./ feed_concentration; 


end
