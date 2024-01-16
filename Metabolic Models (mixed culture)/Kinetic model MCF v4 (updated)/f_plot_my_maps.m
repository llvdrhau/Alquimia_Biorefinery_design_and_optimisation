function [handle] = f_plot_my_maps(compounds,pH,HRT,Qliq_vector,Qgas_vector,states_matrix,inhibition_matrix,AAsugars_comsumptionratio,Vliq,flagPlot)

flag1=flagPlot.flag1;
flag2=logical(flagPlot.flag2);
flag3=logical(flagPlot.flag3);
flag4=logical(flagPlot.flag4);
flag5=flagPlot.flag5;
flag6=logical(flagPlot.flag6);
flag7=logical(flagPlot.flag7);
flag8=logical(flagPlot.flag8);
flag9=logical(flagPlot.flag9);

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
    
fig_number=1;

k=0.4;

%% VFA yield map

if flag2==1
    i=1;
    for i=1:length(pH)
        j=1;
        for j=1:length(HRT)
            VFAyield_matrix(i,j)=100*sum(states_matrix(j+length(HRT)*(i-1),pos_Sac:pos_Sva))/COD_feed;
            j=j+1;
        end
        i=i+1;
    end
    
    for i=2:length(pH)-1
        for j=2:length(HRT)-1
            VFAyield_matrix(i,j)=VFAyield_matrix(i,j)*k+(VFAyield_matrix(i+1,j)+VFAyield_matrix(i-1,j)+VFAyield_matrix(i,j+1)+VFAyield_matrix(i,j-1))*(1-k)/4;
            j=j+1;
        end
        i=i+1;
    end
    
    figure(fig_number)
    f_my_contorno(VFAyield_matrix,pH,HRT)
    title('(b)');
    i=i+1;
else
end


%% Acidification map

if flag3==1
    i=1;
    for i=1:length(pH)
        j=1;
        for j=1:length(HRT)
            aux1=sum(states_matrix(j+length(HRT)*(i-1),pos_Ssu:pos_Sfa));
            aux2=sum(states_matrix(j+length(HRT)*(i-1),pos_Xc1:pos_Xc2));
            aux3=sum(states_matrix(j+length(HRT)*(i-1),pos_Xch:pos_Xi));
            acidification_matrix(i,j)=100*(COD_feed-aux1-aux2-aux3)/COD_feed;
            j=j+1;
        end
        i=i+1;
    end
   
     for i=2:length(pH)-1
        for j=2:length(HRT)-1
            acidification_matrix(i,j)=acidification_matrix(i,j)*k+(acidification_matrix(i+1,j)+acidification_matrix(i-1,j)+acidification_matrix(i,j+1)+acidification_matrix(i,j-1))*(1-k)/4;
            j=j+1;
        end
        i=i+1;
    end
    
   figure (fig_number+1)
    f_my_contorno(acidification_matrix,pH,HRT)
    title('(a)');
    i=i+1;
else
end


%% Mass balance map

if flag4==1
    i=1;
    for i=1:length(pH)
        j=1;
        for j=1:length(HRT)
            Qliq=Qliq_vector(j+length(HRT)*(i-1));
            Qgas=Qgas_vector(j+length(HRT)*(i-1));
            COD_total=sum(states_matrix(j+length(HRT)*(i-1),:));
            notCOD_aux=sum(states_matrix(j+length(HRT)*(i-1),pos_Stic:pos_Sin));
            COD_gas=sum(states_matrix(j+length(HRT)*(i-1),pos_Scat:pos_Sgas_ch4));
            massbalance_matrix(i,j)=100*((COD_total-notCOD_aux-COD_gas)*Qliq+COD_gas*Qgas-COD_feed*Qliq)/COD_feed*Qliq;
            j=j+1;
        end
        i=i+1;
    end
    figure (fig_number+2)
    f_my_contorno(massbalance_matrix,pH,HRT)
    title('Mass balance error (%)');
    i=i+1;
else
end


%
%% VFA spectrum Cimpar/Cpar map

if flag7==1
    pos_Sac=index(strcmp(Abb,'Sac'));
    pos_Sva=index(strcmp(Abb,'Sva'));
    pos_Stic=index(strcmp(Abb,'Stic'));
    pos_Sin=index(strcmp(Abb,'Sin'));
    pos_Scat=index(strcmp(Abb,'Scat'));
    pos_Sgas_ch4=index(strcmp(Abb,'Sgas_ch4'));
    COD_feed=sum(feed)-sum(feed(pos_Stic:pos_Sin))-sum(feed(pos_Scat:pos_Sgas_ch4));
        i=1;
    for i=1:length(pH)
        j=1;
        for j=1:length(HRT)
            Sac=states_matrix(j+length(HRT)*(i-1),pos_Sac)/64;
            Spro=states_matrix(j+length(HRT)*(i-1),pos_Sac+1)/112;
            Sbu=states_matrix(j+length(HRT)*(i-1),pos_Sac+2)/160;
            Sva=states_matrix(j+length(HRT)*(i-1),pos_Sva)/208;
            if Sac>Spro
                HV=Sva+Spro;
                HB=Sbu+(Sac-Spro)/2;
            else
                HV=Sva+Sac;
                HB=Sbu;
            end
            spectrum_matrix(i,j)=HV/(HB+HV);
            j=j+1;
        end
        i=i+1;
    end
    
     for i=2:length(pH)-1
        for j=2:length(HRT)-1
            spectrum_matrix(i,j)=spectrum_matrix(i,j)*k+(spectrum_matrix(i+1,j)+spectrum_matrix(i-1,j)+spectrum_matrix(i,j+1)+spectrum_matrix(i,j-1))*(1-k)/4;
            j=j+1;
        end
        i=i+1;
    end
    
    figure(fig_number+3)
    f_my_contorno(spectrum_matrix,pH,HRT)
    title('(d)');
    i=i+1;
else
end


%% VFA productivity map

if flag6==1
    pos_Sac=index(strcmp(Abb,'Sac'));
    pos_Sva=index(strcmp(Abb,'Sva'));
    pos_Stic=index(strcmp(Abb,'Stic'));
    pos_Sin=index(strcmp(Abb,'Sin'));
    pos_Scat=index(strcmp(Abb,'Scat'));
    pos_Sgas_ch4=index(strcmp(Abb,'Sgas_ch4'));
    COD_feed=sum(feed)-sum(feed(pos_Stic:pos_Sin))-sum(feed(pos_Scat:pos_Sgas_ch4));
    i=1;
    for i=1:length(pH)
        j=1;
        for j=1:length(HRT)
            Qliq=Qliq_vector(j+length(HRT)*(i-1));
            productivity_matrix(i,j)=sum(states_matrix(j+length(HRT)*(i-1),pos_Sac:pos_Sva))*Qliq/Vliq;
            j=j+1;
        end
        i=i+1;
    end
    
    for i=2:length(pH)-1
        for j=2:length(HRT)-1
            productivity_matrix(i,j)=productivity_matrix(i,j)*k+(productivity_matrix(i+1,j)+productivity_matrix(i-1,j)+productivity_matrix(i,j+1)+productivity_matrix(i,j-1))*(1-k)/4;
            j=j+1;
        end
        i=i+1;
    end
    
     figure(fig_number+4)
    f_my_contorno(productivity_matrix,pH,HRT)
    title('(c)');
    i=i+1;
else
end


%% VFAs map

if flag8==1
    i=1;
    for i=1:length(pH)
        j=1;
        for j=1:length(HRT)
            VFAs_matrix(i,j)=sum(states_matrix(j+length(HRT)*(i-1),pos_Sac:pos_Sva));
            j=j+1;
        end
        i=i+1;
    end
    
    for i=2:length(pH)-1
        for j=2:length(HRT)-1
            VFAs_matrix(i,j)=VFAs_matrix(i,j)*k+(VFAs_matrix(i+1,j)+VFAs_matrix(i-1,j)+VFAs_matrix(i,j+1)+VFAs_matrix(i,j-1))*(1-k)/4;
            j=j+1;
        end
        i=i+1;
    end
    
    handle=figure(fig_number+5)
    f_my_contorno(VFAs_matrix,pH,HRT)
    title('VFAs (g/L)');
    i=i+1;
else
end


%% Inhibition maps

if flag5==0
else
    i=1;
    for i=1:length(flag5)
        flag5_pos=flag5(i);
        j=1;
        for j=1:length(pH)
            k=1;
            for k=1:length(HRT)
                flag5_matrix(j,k)=inhibition_matrix(k+length(HRT)*(j-1),flag5_pos);
                k=k+1;
            end
            j=j+1;
        end
        fig_number=i+6;
        figure (fig_number)
        f_my_contorno(flag5_matrix,pH,HRT)
        title('Inhbition');
        i=i+1;
    end
end


%% states map

if flag1==0
else
    i=1;
    for i=1:length(flag1)
        flag1_pos=flag1(i);
        j=1;
        for j=1:length(pH)
            k=1;
            for k=1:length(HRT)
                flag1_matrix(j,k)=states_matrix(k+length(HRT)*(j-1),flag1_pos);
                k=k+1;
            end
            j=j+1;
        end
        
        fig_number=length(flag5)+i+6;
        figure (fig_number)
        f_my_contorno(flag1_matrix,pH,HRT)
        figure_name=Abb(flag1_pos);
        title(figure_name);
        i=i+1;
    end
end

%% Methanisation map

if flag9==1
    i=1;
    for i=1:length(pH)
        j=1;
        for j=1:length(HRT)
            Qliq=Qliq_vector(j+length(HRT)*(i-1));
            Qgas=Qgas_vector(j+length(HRT)*(i-1));
            methanisation_matrix(i,j)=100*states_matrix(j+length(HRT)*(i-1),pos_Sgas_ch4)*Qgas/(COD_feed*Qliq);
            j=j+1;
        end
        i=i+1;
    end
    figure (fig_number+1)
    f_my_contorno(methanisation_matrix,pH,HRT)
    title('Methanisation');
    i=i+1;
else
end

end
