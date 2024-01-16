% This script creates data points for a creation of a sorrogate model. for
% this model the variables are: 
% 1) concentration glucose  
% 2) concentration protein (gelatine or casein)  
% 3) pH 

% Respons variable: acetae, propionate, buterate, valerate, hydrogen and
% Ethanol 

% Lucas Van der Hauwaert. CRETUS Institute. University of Santiago de Compostela. 
% Spain. May 2022. Please contact lucas.vanderhauwaert@usc.es if you intend 
% to use this code.

clear 
close all 
clc 

path = genpath(cd());
addpath(path)


% Select the protein 1 = casein, 2 = gelatine, 3= albumin. 0 = ONLY GLUCOSE
                                                  
proteinID = 0; 
saveName = '60_Data_glucose_pH';
excelName = 'Glucose_PH_60_Data_Points.xlsx';
N = 60;  % Number of data points 
%% fixed variables 

HRT= 3; % Fixed at 3 d
% pHRange= linspace(5,8,5) ;  % choose the HRT interval you want to  
% feed_concentration= 10; % substrates COD concentration kg COD/m^-3 % FIND influence of this 

%% redundant 
Despacho='C:\Users\Reboots\OneDrive - Universidade de Santiago de Compostela\Alquimia\Model_Mateo\Kinetic model MCF v4 (updated)\Stoichiometric model\R Generation';
Casa =   'C:\Users\mateo\OneDrive - Universidade de Santiago de Compostela\Mateo\Modelado\Kinetic model MCF v4\Stoichiometric model\R Generation\';
USC =    'C:\Users\lucas.vanderhauwaert\OneDrive - Universidade de Santiago de Compostela\Alquimia\Model_Mateo\Kinetic model MCF v4 (updated)\Stoichiometric model\R Generation\'; 
%USC = 'C:\Users\lucas\OneDrive - Universidade de Santiago de Compostela\Alquimia\Model_Mateo\Kinetic model MCF v4 (updated)\Stoichiometric model\R Generation\'; 
flagCase.computer=USC;
flagCase.simulation=0;

%% i seems a bit useless (may be change it?)
nr=1;

%% Simulation parameters and options   
tEnd=150;   %days
nIteration=min(5e2,tEnd*(24*6));    % one time node per 10 min
t=linspace(0,tEnd,nIteration);
options = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-12);


%% molar data 
%T = readtable('myMCFModelStructure.xlsx','AA profile');
T = xlsread('myMCFModelStructure.xlsx','AA profile');
AA_MolarMass = T(4,:)'; % g/mol 
glu_MolarMass = 180.156 ; % g/mol 

[T2] = xlsread('Aminoacidos.xlsx','AA'); % Reading AA profile from excel file

if proteinID == 0
    Prot = T2(:,1)' * 0;  % not interested in proteins if ID is 0 
else
    Prot=T2(:,proteinID)'; % Select the protein 1 = casein, 2 = gelatine, 3= albumin.
end 
%%  start loop 
% flagCase structure is refered to the substrate employed in the
% fermentation 
% 1- GELATINE        2- Vinasses       3- Regrind pasta 
% 4- Wheat bran     5- Bread crust    6- Fruit and vegetables waste 
% 7- Tuna canning wastewater

%% select substrate and pH
% one signle substrate as protein input.  so no mixing of protein  sources 
% i.e., flagCase.substrates_ratio = 0;
if proteinID == 0
    flagCase.substrate1 = 1; % arrbatratry:  their is no protein in the system
    flagCase.substrate2 = 1; 
    flagCase.substrates_ratio = 0; % ratio concentration S2/S1 [0:1] 0 if not a mixture
else
    flagCase.substrate1 = proteinID; % choose substrate e.g, GELATINE
    flagCase.substrate2 = proteinID; % choose (optional) 2nd  protein substrate (not used here)
    flagCase.substrates_ratio = 0; % ratio concentration S2/S1 [0:1] 0 if not a mixture
end
%% preallocattion 
features = zeros(N, 3);
Yhat = zeros(N, 7); % 7 outputs 
%% run loop
rng(1)

for i = 1:N
  
    pH = 4 + (9-4).*rand();

    if proteinID == 0
        conc_prot  = 0; % if    o nly glucose is of interest 
    else 
        conc_prot = 1+ (10-1).*rand() ;
    end 

    conc_glu = 10 ;
    
    features(i,:) = [conc_prot, conc_glu, pH];
    
    %% Reading from Excel file
    
    % This function loads the:
    % compounds index, names and initial concentrations
    % feeding concentrations
    % kinetic parameters
    % stoichiometric parameters
    % stoichiometric matrix
    
    [compounds, parameters, stoiMatrix, reaction, AA_profile,...
        new_sol_particulate, new_charb_prot_scarb, new_AA] ...
        = f_loadModelXlsx_MCF(flagCase,nr,10); % 10 doesn't really matter here, we want the constant values here namely mol percentage and MW of AA
    
    %% get the imputs into mol/L 
    % Glucose 
    M_glu = conc_glu / glu_MolarMass ; 
    
    % protein
%     mol_percent_AA = AA_profile(1,:); % mol fractions 
%     MW_selected_protein = mol_percent_AA * AA_MolarMass; % in g/mol 
%     M_Vector_AA = (conc_prot/MW_selected_protein) .* mol_percent_AA;  % in mol/L (M)
    MolarMassAA = [174	89	133	146	147	105	119	89	75	115	117	131	131	117	146	132	155];
    sumProt = sum(Prot);
    percentMM = zeros(length(MolarMassAA),1);

    for j = 1:length(MolarMassAA)
        percentMM(j) = Prot(j)/sumProt;
    end
    MM_Protein = MolarMassAA * percentMM;

    mAA = Prot.*5;
    overallMM_protein = sum(mAA);
    disp(overallMM_protein)
    % in this order 
    % Arg  Ala	Asp	Lys	Glut Ser Thr Cys Gly	Pro	Val	IsoL	Leu	Meth	GluM	AspG	Hist

    M_Vector_AA = Prot/5*conc_prot;
    
    %% Reactor operation parameters
    
    Qliq=1./HRT; %m^3 days^-1
    Vliq = parameters.reactorPar(strcmp(parameters.reactorNames,'Vliq')); % m^3
    
    HRT=Vliq./Qliq;  % days
    D=1./HRT;  % dilution rate days^-1
    
    %% Simulation
    xInitial=compounds.initial;
     
    p=genpath(cd());
    addpath(p)
    load(strcat(flagCase.computer,'R.mat'))    %??
    load(strcat(flagCase.computer,'Opt.mat'))  %??

    % begin with false error lable
    R.errorLable  = false;
    
    Opt.electronSource=0;  % 0 = normal fermentation; 1 = electrofermentation
    Opt.eSource=-0.2;      % offset of the optimiser for NADH balance to reflect the incorporation of wirhdrawal of electron (mol NADH/Lx h)
    
    %% run bioenergetic/kinetic model so stoichiometric table can be updated
    
    [~, End_Conc] = f_determine_my_F(R, Opt, pH, M_Vector_AA, M_glu);
    Yhat(i,:) = End_Conc;  % In mol/L ?? 

%     if R.errorLable == true
%         Yhat(i,:) = End_Conc;
%     else 
%         Yhat(i,:) = End_Conc;
%     end
    
end
%rowCounter = rowCounter + length(pHRange);

%%
nameOutputs = {'acetate_mol/L'
    'propionate_mol/L'
    'butyrate_mol/L'
    'valerate_mol/L'
    'Hydrogen_mol/L'
    'Ethanol_mol/L'
    'Biomass mol/L'
    };

nameInputs = {'concentration_protein_g/L' 'concentration_glucose_g/L' 'pH'}; % in g/L 

%%
save(saveName,'nameOutputs', 'nameInputs','features', 'Yhat');

%% write to Excel 
f_data2excel(features, nameInputs, excelName, 'inputs')
f_data2excel(Yhat, nameOutputs, excelName, 'outputs')

