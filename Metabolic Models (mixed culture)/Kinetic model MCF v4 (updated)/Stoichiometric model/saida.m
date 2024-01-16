function [output,feeding,Cin,C,conversions,substrates] = saida()

horas=0.1;
close all
contador=1;
output=zeros(contador,8);
substrates=zeros(contador,18);
feeding=zeros(contador,17);
for m=1:contador
    name=horzcat('M_',num2str(m),'.mat');
    load(name)
    name2=horzcat('All_StatesVar_',num2str(m),'.mat');
    load(name2)
    %      load('Opt.mat')
    %      load('R.mat')
    R=M.R;
    Opt=M.Opt;
    %%
    datos=All_StatesVar;
    nSt=R.St.numSt;
    y=datos(:,2:nSt+1);
    
    Ac=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ce_Ac')));
    Pro=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ce_Pro')));
    iBut=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ce_iBu')));
    But=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ce_Bu')));
    iVal=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ce_iVal')));
    Val=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ce_Val')));
    iCap=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ce_iCap')));
    For=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ce_For')));
    Et=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ce_EtOH')));
    
    output(m,1)=Ac;
    output(m,2)=Pro;
    output(m,3)=iBut;
    output(m,4)=But;
    output(m,5)=iVal;
    output(m,6)=Val;
    output(m,7)=iCap;
    output(m,8)=Et;
    output(m,9)=For;
    
    Glu=find(strcmp(R.St.StNames,'Ce_Glu')==1);
    AA1=find(strcmp(R.St.StNames,'Ce_Arg')==1);
    AAend=find(strcmp(R.St.StNames,'Ce_Hist')==1);
    pos=[Glu AA1:AAend];
    substrates(m,:)=mean(y(end-horas*10:end,pos));
    
    
    Ac_i=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ci_Ac')));
    Pro_i=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ci_Pro')));
    iBut_i=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ci_iBut')));
    But_i=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ci_Bu')));
    iVal_i=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ci_iVal')));
    Val_i=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ci_Val')));
    iCap_i=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ci_iCap')));
    For_i=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ci_For')));
    Et_i=mean(y(end-horas*10:end,strcmp(R.St.StNames,'Ci_EtOH')));
    
    Cin(m,1)=Ac_i;
    Cin(m,2)=Pro_i;
    Cin(m,3)=iBut_i;
    Cin(m,4)=But_i;
    Cin(m,5)=iVal_i;
    Cin(m,6)=Val_i;
    Cin(m,7)=iCap_i;
    Cin(m,8)=Et_i;
    Cin(m,9)=For_i;
    
    
    
    nAA=length(Opt.index_nz)-1;
    
    C=cell(nAA+2,2);
    vars = {'Pyr','Arg','Ala','Asp','Lys','Glut','Ser','Thr','Cys','Gly','Prol','Vali','IsoL','Leu','Meth','GluM','AspG','Hist'};
    pos_Ana_Prot=strcmp(R.rm.rmNames,'AnabProt');
    pos_1AA = strcmp(horzcat('Ci_',vars{2}),R.St.StNames);
    Ana_Prot = R.rm.r(pos_Ana_Prot)*R.pOp.SRT*R.pOp.Vx/R.pOp.Vr*(-R.rm.stoM(pos_1AA,pos_Ana_Prot));
    
    pos_Ana_Glu=strcmp(R.rm.rmNames,'Anab');
    Ana_Glu = R.rm.r(pos_Ana_Glu)*R.pOp.SRT*R.pOp.Vx/R.pOp.Vr;
    
    C=vars';
    C(2:end+1,:) = C(1:end,:);
    C(1,:) = {''};
    C{2,1}='Glu';
    C{1,2}='Concentración entrada (M)';
    C{1,3}='Concentración saída (M)';
    C{1,4}='Diferencia (M)';
    C{1,5}='Restando Ana (M)';
    C{1,6}='Concentración ext final (M)';
    C{1,7}='Consumo metabólico (base molar)';
    C{1,8}='Consumo catabólico (base molar)';
    for i = 1:length(vars)
        posi=strcmp(horzcat('Ce_',C{i+1,1}),R.St.StNames);
        C{i+1,2}=R.Inf_.InfV(find(posi==1)+2);
        C{i+1,3}=mean(y(end-horas*10:end,posi));
        C{i+1,4}=C{i+1,2}-C{i+1,3};
        if strcmp(horzcat('Ce_',C{i+1,1}),'Ce_Glu')
            C{i+1,5}=C{i+1,4}-Ana_Glu;
        else
            if (C{i+1,4}-Ana_Prot)>0
                C{i+1,5}=C{i+1,4}-Ana_Prot;
            else
                C{i+1,5}=0;
            end
        end
        C{i+1,6}=y(posi);
        C{i+1,7}=C{i+1,4}/C{i+1,2};
        C{i+1,8}=C{i+1,5}/C{i+1,2};
    end
    % GluM, AspG, Hist and Arg produce other AA
%     C{4,5}=C{4,5}+C{3,5};            %Alanine = Alanine + Arginine
%     C{5,5}=C{5,5}+C{18,5};           %Aspartate = Aspartate + Asparagine
%     C{7,5}=C{7,5}+C{17,5}+C{19,5};   %Glutamate = Glutamate + Glutamine + Histidine
    
    conversion_met_prot=sum([C{3:end,4}])/sum([C{3:end,2}]);
    numC=[6;3;4;6;5;3;4;3;2;5;5;6;6;5;5;4;6]';
    conversion_met_C_prot=sum([C{3:end,4}].*numC)/sum([C{3:end,2}].*numC);
    
    conversion_cat_prot=sum([C{3:end,5}])/sum([C{3:end,2}]);
    numC=[6;3;4;6;5;3;4;3;2;5;5;6;6;5;5;4;6]';
    conversion_cat_C_prot=sum([C{3:end,5}].*numC)/sum([C{3:end,2}].*numC);
    
    PM=[174.2 89.0931 133.11 146.19 147.13 105.0926 119.1192 121.16 75.0666 115.13 117.151 131.1729 131.17 149.21 146.15 132.1179 155.1546];
    conversion_met_mass_prot=sum([C{3:end,4}].*PM)/sum([C{3:end,2}].*PM);
    conversion_cat_mass_prot=sum([C{3:end,5}].*PM)/sum([C{3:end,2}].*PM);
     
    conversion_met_glu=sum([C{2,4}])/sum([C{2,2}]);
    conversion_cat_glu=sum([C{2,5}])/sum([C{2,2}]);
    
    conversions(m,1) = conversion_met_mass_prot;
    conversions(m,2) = conversion_cat_mass_prot;
    conversions(m,3) = conversion_met_glu;
    conversions(m,4) = conversion_cat_glu;
    conversions(m,5) = Ana_Prot;
    conversions(m,6) = Ana_Glu;
    conversions(m,7) = y(end,50)*24.6;
    
    pos_glu=strmatch('Ce_Pyr',R.St.StNames)+1;
    feeding(m,1)=M.R.Inf_.InfV(pos_glu);
    pos1=strmatch('Ce_Arg',R.St.StNames)+2;
    posend=strmatch('Ce_Hist',R.St.StNames)+2;
    feeding(m,2:18)=M.R.Inf_.InfV(pos1:posend);
    feeding(m,19)=M.R.Inf_.InfV(2);
    feeding(m,20)=M.R.pOp.SRT;
    %
    %     NADH(m,1)=mean(y(:,strcmp(R.St.StNames,'Ci_NADH')));
    %     NADH(m,2)=mean(y(:,strcmp(R.St.StNames,'Ci_NAD')));
    %     NADH(m,3)=NADH(m,2)/NADH(m,1);
    %     NADH(m,4)=8.314e-3*298.15*abs(log(10)-log(NADH(m,3)));
    %     NADH(m,5)=mean(datos(:,377));
    %     NADH(m,6)=(tanh(-10*(datos(end,160)+2))+1)/2;
    
    %     figure
    %     a=y(:,strcmp(R.St.StNames,'Ci_NAD'))./y(:,strcmp(R.St.StNames,'Ci_NADH'));
    %     plot(datos(:,1),a)
    %     figure
    %     plot(datos(:,1),datos(:,379))
    
    y=datos(:,2:nSt+1);
    Ac=y(:,strcmp(R.St.StNames,'Ce_Ac'));
    Pro=y(:,strcmp(R.St.StNames,'Ce_Pro'));
    But=y(:,strcmp(R.St.StNames,'Ce_Bu'));
    iBut=y(:,strcmp(R.St.StNames,'Ce_iBu'));
    Val=y(:,strcmp(R.St.StNames,'Ce_Val'));
    iVal=y(:,strcmp(R.St.StNames,'Ce_iVal'));
    iCap=y(:,strcmp(R.St.StNames,'Ce_iCap'));
    Et=y(:,strcmp(R.St.StNames,'Ce_EtOH'));
    For=y(:,strcmp(R.St.StNames,'Ce_For'));
    t=datos(:,1);
    figure
    plot(t,Ac,'LineWidth',2)
    hold on
    plot(t,Pro,'LineWidth',2)
    hold on
    plot(t,But,'LineWidth',2)
    hold on
    plot(t,iBut,'LineWidth',2)
    hold on
    plot(t,Val,'LineWidth',2)
    hold on
    plot(t,iVal,'LineWidth',2)
    hold on
    plot(t,iCap,'LineWidth',2)
    hold on
    plot(t,Et,'r','LineWidth',2)
    hold on
    plot(t,For,'k','LineWidth',2)
    hold off
    xlabel('Time(h)','FontSize',12,'FontWeight','bold')
    ylabel('Concentration (M)','FontSize',12,'FontWeight','bold')
    title('Extracellular VFA concentration')
    leg_names={'Acetate','Propionate','Butyrate','iso-Butyrate','Valerate','iso-Valerate','iso-Caproate','Ethanol','Formate'};
    legend(leg_names,'Location','eastoutside','Orientation','vertical','Fontsize',14)
    
%     Ac=y(:,strcmp(R.St.StNames,'Ci_Ac'));
%     Pro=y(:,strcmp(R.St.StNames,'Ci_Pro'));
%     But=y(:,strcmp(R.St.StNames,'Ci_Bu'));
%     iBut=y(:,strcmp(R.St.StNames,'Ci_iBut'));
%     Val=y(:,strcmp(R.St.StNames,'Ci_Val'));
%     iVal=y(:,strcmp(R.St.StNames,'Ci_iVal'));
%     iCap=y(:,strcmp(R.St.StNames,'Ci_iCap'));
%     Et=y(:,strcmp(R.St.StNames,'Ci_EtOH'));
%     For=y(:,strcmp(R.St.StNames,'Ci_For'));
%     t=datos(:,1);
%     figure
%     plot(t,Ac,'LineWidth',2)
%     hold on
%     plot(t,Pro,'LineWidth',2)
%     hold on
%     plot(t,But,'LineWidth',2)
%     hold on
%     plot(t,iBut,'LineWidth',2)
%     hold on
%     plot(t,Val,'LineWidth',2)
%     hold on
%     plot(t,iVal,'LineWidth',2)
%     hold on
%     plot(t,iCap,'LineWidth',2)
%     hold on
%     plot(t,Et,'r','LineWidth',2)
%     hold on
%     plot(t,For,'LineWidth',2)
%     hold off
%     xlabel('Time(h)','FontSize',12,'FontWeight','bold')
%     ylabel('Concentration (M)','FontSize',12,'FontWeight','bold')
%     title('Intracellular VFA concentration')
%     leg_names={'Acetate','Propionate','Butyrate','iso-Butyrate','Valerate','iso-Valerate','iso-Caproate','Ethanol','Formate'};
%     legend(leg_names,'Location','eastoutside','Orientation','vertical','Fontsize',14)
%     
%     figure
%     plot(t,y(:,51))
%     title('External glucose concentration (M)')
%     
%     figure
%     plot(t,datos(:,[216 217]))
%     title('Anabolism and decay rates (mol/L h')
%     legend('Glicosa','Proteína')
end

    %% Datos modificables segundo o caso
    nz=length(Opt.z);
    nSt=R.St.numSt;
    nr=R.rm.num_r;
    nAA=length(Opt.index_nz)-1;
    ruta=Opt.RUTA(1:13,nr-22:nr-6);
    pos_CiH=nSt+2;
    datos=All_StatesVar;
    index=Opt.index_nz;
    zetas=datos(:,end-nz:end-1);
    gibbs=datos(:,nSt+4:nSt+3+nr);
    n_reac=R.rm.rmNames;
    %t_end=size(zetas,1)/60;
    t=datos(:,1);

    %%
    vars = {'Pyr','Arg','Ala','Asp','Lys','Glut','Ser','Thr','Cys','Gly','Prol','Vali','IsoL','Leu','Meth','GluM','AspG','Hist'};

    %% Plotting VFA

    y=datos(:,2:nSt+1);
    Ac=y(:,strcmp(R.St.StNames,'Ce_Ac'));
    Pro=y(:,strcmp(R.St.StNames,'Ce_Pro'));
    But=y(:,strcmp(R.St.StNames,'Ce_Bu'));
    iBut=y(:,strcmp(R.St.StNames,'Ce_iBu'));
    Val=y(:,strcmp(R.St.StNames,'Ce_Val'));
    iVal=y(:,strcmp(R.St.StNames,'Ce_iVal'));
    iCap=y(:,strcmp(R.St.StNames,'Ce_iCap'));
    Et=y(:,strcmp(R.St.StNames,'Ce_EtOH'));
    leg_names={'Acetate','Propionate','Butyrate','iso-Butyrate','Valerate','iso-Valerate','iso-Caproate','Ethanol'};
    F={};
    F(:,1)=leg_names;
    F(2:end+1,:) = F(1:end,:);
    F(1,:) = {''};
    F{1,2}='End OUT Conc (M)';
    F{2,2}=Ac(end);
    F{3,2}=Pro(end);
    F{4,2}=But(end);
    F{5,2}=iBut(end);
    F{6,2}=Val(end);
    F{7,2}=iVal(end);
    F{8,2}=iCap(end);
    F{9,2}=Et(end);

    Ac=y(:,strcmp(R.St.StNames,'Ci_Ac'));
    Pro=y(:,strcmp(R.St.StNames,'Ci_Pro'));
    But=y(:,strcmp(R.St.StNames,'Ci_Bu'));
    iBut=y(:,strcmp(R.St.StNames,'Ci_iBut'));
    Val=y(:,strcmp(R.St.StNames,'Ci_Val'));
    iVal=y(:,strcmp(R.St.StNames,'Ci_iVal'));
    iCap=y(:,strcmp(R.St.StNames,'Ci_iCap'));
    Et=y(:,strcmp(R.St.StNames,'Ci_EtOH'));

    F{1,3}='End IN Conc (M)';
    F{2,3}=Ac(end);
    F{3,3}=Pro(end);
    F{4,3}=But(end);
    F{5,3}=iBut(end);
    F{6,3}=Val(end);
    F{7,3}=iVal(end);
    F{8,3}=iCap(end);
    F{9,3}=Et(end);



    C=cell(nAA+2,2);
    pos_Ana_Pro=strcmp(R.rm.rmNames,'AnabProt');
    pos_BM = strcmp(R.St.StNames,'X');
    coeff_BM=R.rm.stoM(pos_BM,pos_Ana_Pro);
    Ana=mean(y(end-horas*10:end,pos_BM))/coeff_BM/nAA;
    C=vars';
    C(2:end+1,:) = C(1:end,:);
    C(1,:) = {''};
    C{2,1}='Glu';
    C{1,2}='Concentración entrada (M)';
    C{1,3}='Concentración saída (M)';
    C{1,4}='Diferencia (M)';
    C{1,5}='Restando Ana (M)';
    C{1,6}='Concentración ext final (M)';
    for i = 1:length(vars)
        posi=strcmp(strcat('Ce_',C{i+1,1}),R.St.StNames);
        C{i+1,2}=R.mFd.FdV(find(posi==1)+2);
        C{i+1,3}=mean(y(end-horas*10:end,posi));
        C{i+1,4}=C{i+1,2}-C{i+1,3};
        C{i+1,5}=C{i+1,4}-Ana;
        C{i+1,6}=y(end,posi);
    end
    % GluM, AspG, Hist and Arg produce other AA
    C{4,5}=C{4,5}+C{3,5};            %Alanine = Alanine + Arginine
    C{5,5}=C{5,5}+C{18,5};           %Aspartate = Aspartate + Asparagine
    C{7,5}=C{7,5}+C{17,5}+C{19,5};   %Glutamate = Glutamate + Glutamine + Histidine

    conversion=sum([C{3:end,4}])/sum([C{3:end,2}]);
    numC=[6;3;4;6;5;3;4;3;2;5;5;6;6;5;5;4;6]';
    peso=[174;89;133;146;147;105;119;121;75;115;117;131;131;149;146;132;155]';
    conversion_C=sum([C{3:end,4}].*numC)/sum([C{3:end,2}].*numC);
    conversion_M=sum([C{3:end,4}].*peso)/sum([C{3:end,2}].*peso);
    prot_consumida=sum([C{3:end,4}].*peso);
    Yield=mean(All_StatesVar(end-horas*10:end,end));


    output(m,1)=R.pKt.KTr.M_actV(20);
    output(m,2:8)=cell2mat(F(2:8,2));
    output(m,9)=y(end,strcmp(R.St.StNames,'Ce_NH3'));
    for i=1:7
        output(m,i+9)=output(m,i+1)./sum(output(m,2:8));
    end
    output(m,17:23)=cell2mat(F(2:8,3));
    output(m,24)=y(end,strcmp(R.St.StNames,'Ci_NH3'));
    output(m,25)=conversion;
    output(m,26)=conversion_C;
    output(m,27)=conversion_M;
    output(m,28)=Yield;
    pos_ATP=strcmp(R.St.StNames,'Ci_ATP');
    output(m,29)=sum(R.rm.stoM(pos_ATP,6:end).*R.rm.r(6:end)');
    output(m,30)=output(m,27)*R.pOp.Vx;
    output(m,31)=R.pOp.Xt*24.6;
% end
end