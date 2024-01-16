% This script launchs the simulation of the production of VFAs from complex
% organic matter in a CSTR.

% The model is based on AMD1 and has 32 different states and 21
% reactions. For the estimation of the acidogenic fermentation
% stoichiometric parameters, a bioenergetic model developed by Alberte 
% Regueira López is employed.

% Mateo Saavedra del Oso. CRETUS Institute. University of Santiago de Compostela. 
% Spain. November 2020. Please contact msaavedra.deloso@usc.es if you intend 
% to use this code.



%% Options menu and instructions
 
HRT_interval=linspace(2,6,5); % choose the HRT interval you want to simulate
pH_interval=linspace(5,8,5);  % choose the pH interval you want to  
feed_concentration=77.4; % substrates COD concentration kg COD/m^-3

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
flagCase.substrate1=15; % choose substrate 1
flagCase.substrate2=7; % choose substrate 2
flagCase.substrates_ratio=0; % ratio concentration S2/S1 [0:1]
flagCase.simulation=1;
Despacho='C:\Users\Reboots\OneDrive - Universidade de Santiago de Compostela\Alquimia\Model_Mateo\Kinetic model MCF v4 (updated)\Stoichiometric model\R Generation\';
%Casa='C:\Users\mateo\OneDrive - Universidade de Santiago de Compostela\Mateo\Modelado\Kinetic model MCF v4 (updated)\Stoichiometric model\R Generation\';
flagCase.computer=Despacho;
i=1;



% flag plot structure
% -flag1: choose which states would you like to represent, please
% introduce the state or states vector index (see below)
% 1- 'Ssu'   2- 'Saa'    3- 'Sfa'    4- 'Sac'    5- 'Spro'   6- 'Sbu'        
% 7- 'Sva'   8- 'Sh2'    9- 'Set'    10- 'Sch4'  11-'Stic'   12- 'Sin'
% 13- 'Ssi'  14- 'Xc1'   15- 'Xc2'   16- 'Xc3'   17- 'Xch'   18- 'Xpr'
% 19- 'Xtdf' 20- 'Xli'   21- 'Xi'    22- 'Xaf'   23- 'Xfa'   24- 'Xc4'
% 25- 'Xpro' 26- 'Xac'   27- 'Xh2'   28- 'Scat'  29- 'San'   30- 'Sgas_h2'
% 31- 'Sgas_co2'         32- 'Sgas_ch4'

% -flag2: activate the VFA yield plot (1) or not (0)

% -flag3: activate the acidification plot (1) or not (0)

% -flag4: activate the mass balance error plot (1) or not (0)

% -flag5: choose which inhibition would you like to represent (see
% below)
% 1- Iin_lim 2- Ih2_fa 3- Ih2_c4 4- Ih2_pro 5- Inh3 6- IpH_aa 
% 7- IpH_ac  8- IpH_h2 

% -flag6: activate the VFA productivity plot (1) or not (0)

% -flag7: activate the ratio VFA C impar/C par plot (1) or not (0)

flagPlot.flag1=0;
flagPlot.flag2=0;
flagPlot.flag3=0;
flagPlot.flag4=0;
flagPlot.flag5=0;
flagPlot.flag6=1;
flagPlot.flag7=0;
flagPlot.flag8=0;
flagPlot.flag9=0;
%% Reading from Excel file

% This function loads the:
% compounds index, names and initial concentrations
% feeding concentrations
% kinetic parameters
% stoichiometric parameters
% stoichiometric matrix

[compounds, parameters, stoiMatrix, reaction, AA_profile] = f_loadModelXlsx_MCF(flagCase,i,feed_concentration);



%% Reactor operation parameters


Qliq=1./HRT_interval; %m^3 days^-1
Vliq = parameters.reactorPar(strcmp(parameters.reactorNames,'Vliq')); % m^3
pH=pH_interval;

HRT=Vliq./Qliq;  % days
D=1./HRT;  % dilution rate days^-1

%% Simulation parameters and options

xInitial=compounds.initial;     
tEnd=150;   %days
nIteration=min(5e2,tEnd*(24*6));    % one time node per 10 min
t=linspace(0,tEnd,nIteration);
options = odeset('NonNegative',1,'RelTol',1e-10,'AbsTol',1e-12);


%% Simulation

if flagCase.simulation==1
    pH_vector=zeros(length(D)*length(pH),1);
    Qliq_vector=zeros(length(D)*length(pH),1);
    Qgas_vector=zeros(length(D)*length(pH),1);
    states_matrix=zeros(length(D)*length(pH),length(compounds.Abb));
    inhibition_matrix=zeros(length(D)*length(pH),8);
    F_matrix=zeros(length(D)*length(pH),6);
    tsimulation_vector=zeros(length(pH),1);
    AAsugars_comsumptionratio=zeros(length(D)*length(pH),1);
    F_comparison=zeros(length(pH)*2,6);
    i=1;
    for i=1:length(pH)
        tic
        j=1;
        Df=D(round(length(D)/2));
        pHi=pH(i);
        [t,states]=ode15s(@f_my_reactor,t,xInitial,options,compounds,parameters, stoiMatrix,pHi,Df);
        [consumed_AA, consumed_glucose] = f_determine_my_consumed_AA_glucose(states, compounds, parameters,AA_profile,flagCase);
        p=genpath(cd());
         addpath(p)
         load(strcat(flagCase.computer,'R.mat'))
         load(strcat(flagCase.computer,'Opt.mat'))
         Opt.electronSource=0;  % 0 = normal fermentation; 1 = electrofermentation
         Opt.eSource=-0.2;      % offset of the optimiser for NADH balance to reflect the incorporation of wirhdrawal of electron (mol NADH/Lx h)
         [F_acidogenic_fermentation] = f_determine_my_F(R, Opt, pHi, consumed_AA, consumed_glucose);
         Yaf = parameters.stoiPar(strcmp(parameters.stoiNames,'Yaf'));             % yield acidogenic fermentation biomass
         AcidificationSSu_index=cell2mat(reaction.Index(strcmp(reaction.Names,'Acidogenic fermentation Su')));
         AcidificationSaa_index=cell2mat(reaction.Index(strcmp(reaction.Names,'Acidogenic fermentation AA')));
         Sac_index=compounds.index(strcmp(compounds.Abb,'Sac'));
         Set_index=compounds.index(strcmp(compounds.Abb,'Set'));
         stoiMatrix(Sac_index:Set_index,AcidificationSSu_index)=(1-Yaf)*F_acidogenic_fermentation;
         stoiMatrix(Sac_index:Set_index,AcidificationSaa_index)=(1-Yaf)*F_acidogenic_fermentation;
        for j=1:length(D)
            Di=D(j);
            [t,states]=ode15s(@f_my_reactor,t,xInitial,options,compounds,parameters, stoiMatrix,pHi,Di);
            [consumed_AA, consumed_glucose] = f_determine_my_consumed_AA_glucose(states, compounds, parameters,AA_profile,flagCase);
            for k=1:length(t)
                y=states(k,:)';
                Qgas(k)=f_my_gasPhase(y,compounds, parameters);
            end
            for h=1:length(t)
                Sin_index=compounds.index(strcmp(compounds.Abb,'Sin'));
                Sh2_index=compounds.index(strcmp(compounds.Abb,'Sh2'));
                Sin=states(h,Sin_index);
                Sh2=states(h,Sh2_index);
                Inhibition(h,:)=f_my_inhibitions(Sin,Sh2,parameters, pHi);
            end
            pH_vector(j+length(D)*(i-1),1)=pHi;
            Qliq_vector(j+length(D)*(i-1),1)=Vliq*D(j);
            Qgas_vector(j+length(D)*(i-1))=Qgas(end);
            states_matrix(j+length(D)*(i-1),:)= states(end,:);
            inhibition_matrix(j+length(D)*(i-1),:)=Inhibition(end,:);
            F_matrix(j+length(D)*(i-1),:)=F_acidogenic_fermentation';
            AAsugars_comsumptionratio(j+length(D)*(i-1))=sum(consumed_AA)/consumed_glucose;
            j=j+1;
        end
        toc
        tsimulation_vector(i)=toc;
        i=i+1;
    end
    
elseif flagCase.simulation==2
    pH_vector=zeros(length(flagCase.substrates_ratio)*length(pH),1);
    Qliq_vector=zeros(length(flagCase.substrates_ratio)*length(pH),1);
    Qgas_vector=zeros(length(flagCase.substrates_ratio)*length(pH),1);
    states_matrix=zeros(length(flagCase.substrates_ratio)*length(pH),length(compounds.Abb));
    inhibition_matrix=zeros(length(flagCase.substrates_ratio)*length(pH),8);
    F_matrix=zeros(length(flagCase.substrates_ratio)*length(pH),6);
    substrates_mixture_vector=zeros(length(flagCase.substrates_ratio)*length(pH),1);
    i=1;
    for i=1:length(flagCase.substrates_ratio)
        j=1;
        for j=1:length(pH)
            pHj=pH(j);
            flagCase.substrates_ratioj=flagCase.substrates_ratio(i);
            [compounds, parameters, stoiMatrix, reaction, AA_profile] = f_loadModelXlsx_MCF(flagCase, i, feed_concentration);
            [t,states]=ode15s(@f_my_reactor,t,xInitial,options,compounds,parameters, stoiMatrix,pHj,D);
            [consumed_AA, consumed_glucose] = f_determine_my_consumed_AA_glucose(states, compounds, parameters,AA_profile,flagCase);
            p=genpath(cd());
            addpath(p)
            load(strcat(flagCase.computer,'R.mat'))
            load(strcat(flagCase.computer,'Opt.mat'))
            Opt.electronSource=0;  % 0 = normal fermentation; 1 = electrofermentation
            Opt.eSource=-0.2;      % offset of the optimiser for NADH balance to reflect the incorporation of wirhdrawal of electron (mol NADH/Lx h)
             [F_acidogenic_fermentation] = f_determine_my_F(R, Opt, pHj, consumed_AA, consumed_glucose);
             Yaf = parameters.stoiPar(strcmp(parameters.stoiNames,'Yaf'));             % yield acidogenic fermentation biomass
             AcidificationSSu_index=cell2mat(reaction.Index(strcmp(reaction.Names,'Acidogenic fermentation Su')));
             AcidificationSaa_index=cell2mat(reaction.Index(strcmp(reaction.Names,'Acidogenic fermentation AA')));
             Sac_index=compounds.index(strcmp(compounds.Abb,'Sac'));
             Set_index=compounds.index(strcmp(compounds.Abb,'Set'));
             stoiMatrix(Sac_index:Set_index,AcidificationSSu_index)=(1-Yaf)*F_acidogenic_fermentation;
             stoiMatrix(Sac_index:Set_index,AcidificationSaa_index)=(1-Yaf)*F_acidogenic_fermentation;
             [t,states]=ode15s(@f_my_reactor,t,xInitial,options,compounds,parameters, stoiMatrix,pHj,D);
            for k=1:length(t)
                y=states(k,:)';
                Qgas(k)=f_my_gasPhase(y,compounds, parameters);
            end
            for h=1:length(t)
                Sin_index=compounds.index(strcmp(compounds.Abb,'Sin'));
                Sh2_index=compounds.index(strcmp(compounds.Abb,'Sh2'));
                Sin=states(h,Sin_index);
                Sh2=states(h,Sh2_index);
                Inhibition(h,:)=f_my_inhibitions(Sin,Sh2,parameters, pHj);
            end
            pH_vector(j+length(pH)*(i-1),1)=pHj;
            substrates_mixture_vector(j+length(pH)*(i-1),1)=flagCase.substrates_ratioj;
            Qliq_vector(j+length(pH)*(i-1),1)=Vliq*D;
            Qgas_vector(j+length(pH)*(i-1))=Qgas(end);
            states_matrix(j+length(pH)*(i-1),:)= states(end,:); % of interset
            inhibition_matrix(j+length(pH)*(i-1),:)=Inhibition(end,:);
            F_matrix(j+length(pH)*(i-1),:)=F_acidogenic_fermentation';
            j=j+1;
        end
        i=i+1;
    end
end


%% Print the data

save('tunapasta_50x50')

%!shutdown /p
% f_plot_my_maps(compounds,pH,HRT,Qliq_vector,Qgas_vector,states_matrix,inhibition_matrix,AAsugars_comsumptionratio,Vliq,flagPlot)