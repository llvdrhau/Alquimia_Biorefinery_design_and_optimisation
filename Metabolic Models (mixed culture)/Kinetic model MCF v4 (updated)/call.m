% This script launchs the simulation of the production of VFAs from complex
% organic matter in a CSTR.

% The model is based on AMD1 and has 32 and 21 different states and
% reactions respectively. For the estimation of the acidogenic fermentation
% stoichiometric parameters, a bioenergetic model developed by Alberte 
% Regueira López is employed.

% Mateo Saavedra del Oso. CRETUS Institute. University of Santiago de Compostela. 
% Spain. November 2020. Please contact msaavedra.deloso@usc.es if you intend 
% to use this code.



%% Options menu and instructions
 
HRT=4; % choose the HRT interval you want to simulate
pH=7;  % choose the HRT interval you want to  
feed_concentration=10; % substrates COD concentration kg COD/m^-3

% flagCase structure is refered to the substrate employed in the
% fermentation 
% 1- Spinach        2- Vinasses       3- Regrind pasta 
% 4- Wheat bran     5- Bread crust    6- Fruit and vegetables waste 
% 7- Tuna canning wastewater

% and the simulation options: 
% -Option 'Quick' only runs 1 time the bioenergetic model per pH 
% -Option 'Substrates mixture' when 2 substrates are employed at a fixed HRT 
% and the bioenergetic model runs 1 time per pH
% -Else only 1 substrate is employed and bioenergetic model runs per each
% HRT and pH
flagCase.substrate1=7; % choose substrate 1
flagCase.substrate2=7; % choose substrate 2
flagCase.substrates_ratio=0; % ratio concentration S2/S1 [0:1]
Despacho='C:\Users\Reboots\OneDrive - Universidade de Santiago de Compostela\Alquimia\Model_Mateo\Kinetic model MCF v4 (updated)\Stoichiometric model\R Generation\';

Casa='C:\Users\mateo\OneDrive - Universidade de Santiago de Compostela\Mateo\Modelado\Kinetic model MCF v4\Stoichiometric model\R Generation\';

flagCase.computer=Despacho;
flagCase.simulation=0;

i=1;




%% Reading from Excel file

% This function loads the:
% compounds index, names and initial concentrations
% feeding concentrations
% kinetic parameters
% stoichiometric parameters
% stoichiometric matrix

[compounds, parameters, stoiMatrix, reaction, AA_profile] = f_loadModelXlsx_MCF(flagCase,i,feed_concentration);



%% Reactor operation parameters


Qliq=1./HRT; %m^3 days^-1
Vliq = parameters.reactorPar(strcmp(parameters.reactorNames,'Vliq')); % m^3

HRT=Vliq./Qliq;  % days
D=1./HRT;  % dilution rate days^-1

%% Simulation parameters and options

xInitial=compounds.initial;     
tEnd=150;   %days
nIteration=min(5e2,tEnd*(24*6));    % one time node per 10 min
t=linspace(0,tEnd,nIteration);
options = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-12);


%% Simulation

[t,states]=ode15s(@f_my_reactor,t,xInitial,options,compounds,parameters, stoiMatrix,pH,D);
[consumed_AA, consumed_glucose] = f_determine_my_consumed_AA_glucose(states, compounds, parameters,AA_profile, flagCase);

p=genpath(cd());
addpath(p)
load(strcat(flagCase.computer,'R.mat'))
load(strcat(flagCase.computer,'Opt.mat'))
Opt.electronSource=0;  % 0 = normal fermentation; 1 = electrofermentation
Opt.eSource=-0.2;      % offset of the optimiser for NADH balance to reflect the incorporation of wirhdrawal of electron (mol NADH/Lx h)
%%
[F_acidogenic_fermentation] = f_determine_my_F(R, Opt, pH, consumed_AA, consumed_glucose);

%%
Yaf = parameters.stoiPar(strcmp(parameters.stoiNames,'Yaf'));             % yield acidogenic fermentation biomass
AcidificationSSu_index=cell2mat(reaction.Index(strcmp(reaction.Names,'Acidogenic fermentation Su')));
AcidificationSaa_index=cell2mat(reaction.Index(strcmp(reaction.Names,'Acidogenic fermentation AA')));
Sac_index=compounds.index(strcmp(compounds.Abb,'Sac'));
Set_index=compounds.index(strcmp(compounds.Abb,'Set'));

stoiMatrix(Sac_index:Set_index,AcidificationSSu_index)=(1-Yaf)*F_acidogenic_fermentation;
stoiMatrix(Sac_index:Set_index,AcidificationSaa_index)=(1-Yaf)*F_acidogenic_fermentation;

tic
[t,states]=ode15s(@f_my_reactor,t,xInitial,options,compounds,parameters, stoiMatrix,pH,D);
toc


%%
f_plot_my_results(t,states, compounds)

%%
for i=1:length(t)
    y=states(i,:)';
    Qgas(i)=f_my_gasPhase(y,compounds, parameters);
end

for i=1:length(t)
    Sin_index=compounds.index(strcmp(compounds.Abb,'Sin'));
    Sh2_index=compounds.index(strcmp(compounds.Abb,'Sh2'));
    Sin=states(i,Sin_index);
    Sh2=states(i,Sh2_index);
    Inhibition(i,:)=f_my_inhibitions(Sin,Sh2,parameters, pH);
end


figure
plot(t,Qgas)

figure
plot(t,Inhibition)


index=compounds.index;
Abb=compounds.Abb;
feed=compounds.feed;

pos_Ssu=index(strcmp(Abb,'Ssu'));
pos_Sfa=index(strcmp(Abb,'Sfa'));
pos_Xc1=index(strcmp(Abb,'Xc1'));
pos_Xc2=index(strcmp(Abb,'Xc2'));
pos_Xch=index(strcmp(Abb,'Xch'));
pos_Xi=index(strcmp(Abb,'Xi'));
pos_Sac=index(strcmp(Abb,'Sac'));
pos_Sva=index(strcmp(Abb,'Sva'));
pos_Stic=index(strcmp(Abb,'Stic'));
pos_Sin=index(strcmp(Abb,'Sin'));
pos_Scat=index(strcmp(Abb,'Scat'));
pos_Sgas_ch4=index(strcmp(Abb,'Sgas_ch4'));

COD_feed=sum(feed)-sum(feed(pos_Stic:pos_Sin))-sum(feed(pos_Scat:pos_Sgas_ch4));

aux1=sum(states(end,pos_Ssu:pos_Sfa));
aux2=sum(states(end,pos_Xc1:pos_Xc2));
aux3=sum(states(end,pos_Xch:pos_Xi));

acidification = 100*(COD_feed-aux1-aux2-aux3)/COD_feed

VFAyield=100*sum(states(end,pos_Sac:pos_Sva))/COD_feed



