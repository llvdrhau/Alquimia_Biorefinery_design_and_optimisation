
horas=1;
close all
name='1.mat';
name2=horzcat('All_StatesVar_',name);
load(name2)
name3=horzcat('M_',name);
load(name3);
R=M.R;
Opt=M.Opt;
%% Datos modificables segundo o caso
nz=length(Opt.z);
nSt=R.St.numSt;
nr=R.rm.num_r;
nAA=length(Opt.index_nz)-1;
ruta=Opt.RUTA(1:14,nr-23:nr-6);
pos_CiH=nSt+2;
datos=All_StatesVar;
index=Opt.index_nz;
zetas=datos(:,end-nz:end-1);
gibbs=datos(:,nSt+4:nSt+3+nr);
n_reac=R.rm.rmNames;
%t_end=size(zetas,1)/60;
t=datos(:,1);
r=datos(:,nSt+nr+4:nSt+nr+4+nr-1);

%%
vars = {'Pyr','Arg','Ala','Asp','Lys','Glut','Ser','Thr','Cys','Gly','Prol','Vali','IsoL','Leu','Meth','GluM','AspG','Hist'};
Z = [];
Z.(vars{1}) = zetas(:,1:size(ruta,1));
for i = 1:length(vars)-1
  Z.(vars{i+1}) = zetas(:,index(i)+1:index(i+1));
end

G = [];
for i = 1:length(vars)
  G.(vars{i}) =gibbs(:,strmatch(vars{i},n_reac));
end
G.BM=gibbs(:,1:5);
G.Aux=gibbs(:,nr-5:end);
G.Pyr=G.Pyr*ruta';
G.Gibbs=(gibbs(end,:)*Opt.RUTA')';

T=vars';
T(2:end+1,:) = T(1:end,:);
T(1,:) = {''};
T{1,2}='Opcion fav';
T{1,3}='Valor 10h';
T{1,4}='deltaG';

for i=1:length(vars)
a=eval(strcat('Z.',vars{i}));
b=mean(a(end-horas*10:end,:));
[c,d]=max(b);
T{i+1,2}=d;
T{i+1,3}=c;
g=eval(strcat('G.',vars{i}));
h=mean(g(end-horas*10:end,d));
T{i+1,4}=h;

figure
subplot(1,2,1)
plot(t,a,'LineWidth',1.5)
title(strcat('zetas ',vars{i}))
xlabel('Time (h)','FontSize',12,'FontWeight','bold')
ylabel('Zetas','FontSize',12,'FontWeight','bold')
xlim([0 t(end)])
% set(gca,'LineWidth',5,'FontSize',12,'FontWeight','bold') 

subplot(1,2,2)
plot(t,g,'LineWidth',1.5)
title(strcat('\DeltaG ',vars{i}))
xlabel('Time (h)','FontSize',12,'FontWeight','bold')
ylabel('\DeltaG (kJ/mol)','FontSize',12,'FontWeight','bold') 
xlim([0 t(end)])
leg=1:size(a,2);
leg=num2str(leg);
leg=strsplit(leg);
legend(leg,'Location','eastoutside','Orientation','vertical','Fontsize',14)
end
disp(T)

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
plot(t,Et,'LineWidth',2)
hold off
xlabel('Time(h)','FontSize',12,'FontWeight','bold')
ylabel('Concentration (M)','FontSize',12,'FontWeight','bold')
title('VFA concentration')
leg_names={'Acetate','Propionate','Butyrate','iso-Butyrate','Valerate','iso-Valerate','iso-Caproate','Ethanol'};
legend(leg_names,'Location','eastoutside','Orientation','vertical','Fontsize',14)

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

disp(F)

figure
plot(t,-log10(datos(:,pos_CiH)),'LineWidth',2)
xlabel('Time(h)','FontSize',12,'FontWeight','bold')
ylabel('pH','FontSize',12,'FontWeight','bold')
title('Internal pH')

C=cell(nAA+2,2);
    vars = {'Pyr','Arg','Ala','Asp','Lys','Glut','Ser','Thr','Cys','Gly','Prol','Vali','IsoL','Leu','Meth','GluM','AspG','Hist'};
    pos_Ana_Prot=strcmp(R.rm.rmNames,'AnabProt');
    pos_1AA = strcmp(horzcat('Ci_',vars{2}),R.St.StNames);
    Ana_Prot = mean(r(end-horas*10:end,pos_Ana_Prot))*R.pOp.SRT*R.pOp.Vx/R.pOp.Vr*(-R.rm.stoM(pos_1AA,pos_Ana_Prot));
    
    pos_Ana_Glu=strcmp(R.rm.rmNames,'Anab');
    Ana_Glu = mean(r(end-horas*10:end,pos_Ana_Glu))*R.pOp.SRT*R.pOp.Vx/R.pOp.Vr;
    
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

conversion_met_prot=sum([C{3:end,4}])/sum([C{3:end,2}]);
numC=[6;3;4;6;5;3;4;3;2;5;5;6;6;5;5;4;6]';
conversion_met_C_prot=sum([C{3:end,4}].*numC)/sum([C{3:end,2}].*numC);

conversion_cat_prot=sum([C{3:end,5}])/sum([C{3:end,2}]);
numC=[6;3;4;6;5;3;4;3;2;5;5;6;6;5;5;4;6]';
conversion_cat_C_prot=sum([C{3:end,5}].*numC)/sum([C{3:end,2}].*numC);

conversion_met_glu=sum([C{2,4}])/sum([C{2,2}]);
conversion_cat_glu=sum([C{2,5}])/sum([C{2,2}]);

COD = [176 96 96 224 144 80 128 144 48 176 192 240 240 240 144 96 160];
conversion_COD = sum([C{3:end,4}].*COD)/sum([C{3:end,2}].*COD);

fprintf('Protein conversion (C-mol) is: %f \n',conversion_met_C_prot)
fprintf('Glucose conversion (C-mol) is: %f \n',conversion_met_glu)


fprintf('Conversion (COD basis) is: %f \n',conversion_COD)

fprintf('Glucose anabolism (C-molX/C-mol glucose) is: %f \n',Ana_Glu/(6*C{2,4}))
fprintf('Protein anabolism (C-molX/C-mol protein) is: %f \n',full(Ana_Prot/mean(numC'.*(cell2mat(C(3:19,4))))))
fprintf('Protein anabolism (g COD BM/g COD protein) is: %f \n',full(Ana_Prot*1.37*24.6/mean(COD'.*(cell2mat(C(3:19,4))))))

