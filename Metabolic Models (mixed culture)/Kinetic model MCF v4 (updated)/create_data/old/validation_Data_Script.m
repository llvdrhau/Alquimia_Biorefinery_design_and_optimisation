% create validation data 

% This script creates data points for a creation of a sorrogate model. for
% this model the variables are: 
% 1) Fractions of the input (soluble/inert)
% 2) Dissintegration constant kDis
% 3) pH 
% 3) Amino acid compostion 

% respons variable: yields obtained from the bioenergectic model 

% Lucas Van der Hauwaert. CRETUS Institute. University of Santiago de Compostela. 
% Spain. May 2022. Please contact lucas.vanderhauwaert@usc.es if you intend 
% to use this code.

clear 
close all 
clc 

%% fixed variables 

HRT= 3; % Fixed at 3 d
pHRange= linspace(5,8,5) ;  % choose the HRT interval you want to  
feed_concentration= 10; % substrates COD concentration kg COD/m^-3 % FIND influence of this 
N = 20; % number of simulations to run 


%% redundant 
Despacho='C:\Users\Reboots\OneDrive - Universidade de Santiago de Compostela\Alquimia\Model_Mateo\Kinetic model MCF v4 (updated)\Stoichiometric model\R Generation\';
Casa='C:\Users\mateo\OneDrive - Universidade de Santiago de Compostela\Mateo\Modelado\Kinetic model MCF v4\Stoichiometric model\R Generation\';

flagCase.computer=Despacho;
flagCase.simulation=0;

%% i seems a bit useless (may be change it?)
nr=1;

%% Simulation parameters and options   
tEnd=150;   %days
nIteration=min(5e2,tEnd*(24*6));    % one time node per 10 min
t=linspace(0,tEnd,nIteration);
options = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-12);

%%  start loop 

% flagCase structure is refered to the substrate employed in the
% fermentation 
% 1- Spinach        2- Vinasses       3- Regrind pasta 
% 4- Wheat bran     5- Bread crust    6- Fruit and vegetables waste 
% 7- Tuna canning wastewater

combineSubtrateFlags = nchoosek(1:7,2); % make comibination of substrates 

rows = N;  % the number of simulations 
colInputMatrix = 23;
colOutputMatrix = 8;  

inputMatrix = zeros(rows, colInputMatrix); 
outputMatrix = zeros(rows, colOutputMatrix); 
rowCounter = 0; 

rng(1)
for i = 1:N
    pHRange = 5 + (8-5).*rand();
    splitRatio = rand(); %  randomise split
    flagCase.substrates_ratio = splitRatio; % (change for validation set)
    
    %% select substrate and pH
    flagCase.substrate1 = combineSubtrateFlags(i,1); % choose substrate 1
    flagCase.substrate2 = combineSubtrateFlags(i,2); % choose substrate 2
    pH = pHRange;
    %% Reading from Excel file
    
    % This function loads the:
    % compounds index, names and initial concentrations
    % feeding concentrations
    % kinetic parameters
    % stoichiometric parameters
    % stoichiometric matrix
    
    [compounds, parameters, stoiMatrix, reaction, AA_profile,...
        new_sol_particulate, new_charb_prot_scarb, new_AA] ...
        = f_loadModelXlsx_MCF(flagCase,nr,feed_concentration);
    
    %% input characterisation
    inputs = [new_AA, new_sol_particulate, new_charb_prot_scarb, pH ];
    inputMatrix(i,:) = inputs;
    
    %% Reactor operation parameters
    
    Qliq=1./HRT; %m^3 days^-1
    Vliq = parameters.reactorPar(strcmp(parameters.reactorNames,'Vliq')); % m^3
    
    HRT=Vliq./Qliq;  % days
    D=1./HRT;  % dilution rate days^-1
    
    %% Simulation
    xInitial=compounds.initial;
    
    [t,states]=ode15s(@f_my_reactor,t,xInitial,options,compounds,parameters, stoiMatrix,pH,D);
    [consumed_AA, consumed_glucose] = f_determine_my_consumed_AA_glucose(states, compounds, parameters,AA_profile, flagCase);
    
    p=genpath(cd());
    addpath(p)
    load(strcat(flagCase.computer,'R.mat'))    %??
    load(strcat(flagCase.computer,'Opt.mat'))  %??
    
    Opt.electronSource=0;  % 0 = normal fermentation; 1 = electrofermentation
    Opt.eSource=-0.2;      % offset of the optimiser for NADH balance to reflect the incorporation of wirhdrawal of electron (mol NADH/Lx h)
    
    %% run bioenergetic model so stoichiometric table can be updated
    [F_acidogenic_fermentation] = f_determine_my_F(R, Opt, pH, consumed_AA, consumed_glucose);
    
    %% update stoichiometric model
    Yaf = parameters.stoiPar(strcmp(parameters.stoiNames,'Yaf'));             % yield acidogenic fermentation biomass
    AcidificationSSu_index=cell2mat(reaction.Index(strcmp(reaction.Names,'Acidogenic fermentation Su')));
    AcidificationSaa_index=cell2mat(reaction.Index(strcmp(reaction.Names,'Acidogenic fermentation AA')));
    Sac_index=compounds.index(strcmp(compounds.Abb,'Sac'));
    Set_index=compounds.index(strcmp(compounds.Abb,'Set'));
    
    stoiMatrix(Sac_index:Set_index,AcidificationSSu_index)=(1-Yaf)*F_acidogenic_fermentation;
    stoiMatrix(Sac_index:Set_index,AcidificationSaa_index)=(1-Yaf)*F_acidogenic_fermentation;
    
    %% run stoichiometric model and calculate the final yields using
    %  updated info of the stoichiometry
    [t,states]=ode15s(@f_my_reactor,t,xInitial,options,compounds,parameters, stoiMatrix,pH,D);
    
    index=compounds.index;
    Abb=compounds.Abb;
    feed=compounds.feed;
    
    
    pos_Sac=index(strcmp(Abb,'Sac'));
    pos_Sva=index(strcmp(Abb,'Sva'));
    pos_Sch4 = index(strcmp(Abb,'Sch4'));
    
    pos_Stic=index(strcmp(Abb,'Stic'));
    pos_Sin=index(strcmp(Abb,'Sin'));
    pos_Scat=index(strcmp(Abb,'Scat'));
    pos_Sgas_ch4=index(strcmp(Abb,'Sgas_ch4'));
    
    
    
    COD_feed=sum(feed)-sum(feed(pos_Stic:pos_Sin))-sum(feed(pos_Scat:pos_Sgas_ch4));
    
    OutputYields = 100*(states(end,pos_Sac:pos_Sch4))/COD_feed;
    yieldCH4 = 100*(states(end,pos_Sch4))/COD_feed;
    
    allOutputYield = [OutputYields,yieldCH4];
    outputMatrix(i,:) = allOutputYield;
    
end
%%
nameOutputs = {'acetate'
    'propionate'
    'butyrate'
    'valerate'
    'Hydrogen'
    'Ethanol'
    'Methane'
    'mathane_Gas'};

nameInputs = {'Arg'
    'Ala'
    'Asp'
    'Lys'
    'Glut'
    'Ser'
    'Thr'
    'Cys'
    'Gly'
    'Pro'
    'Val'
    'IsoL'
    'Leu'
    'Meth'
    'GluM'
    'AspG'
    'Hist'
    'percentage Particulate matter'
    'percentage soluble'
    'percentage Carbohydrates'
    'percentage protein'
    'percentage short charbohydrates'
    'pH'};

%%
save('Validation_Data_4_ML_model','nameOutputs', 'nameInputs','outputMatrix', 'inputMatrix');
