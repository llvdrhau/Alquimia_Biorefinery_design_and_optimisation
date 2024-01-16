function handle = f_plot_my_maps2(compounds,pH,Qliq_vector,states_matrix,inhibition_matrix,Vliq,flagPlot,flagCase)

flag1=flagPlot.flag1;
flag2=logical(flagPlot.flag2);
flag3=logical(flagPlot.flag3);
flag4=logical(flagPlot.flag4);
flag5=flagPlot.flag5;
flag6=logical(flagPlot.flag6);
flag7=logical(flagPlot.flag7);
flag8=logical(flagPlot.flag8);

substrates_ratio=flagCase.substrates_ratio;

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
    for i=1:length(substrates_ratio)
        j=1;
        for j=1:length(pH)
            VFAyield_matrix(j,i)=100*sum(states_matrix(j+length(pH)*(i-1),pos_Sac:pos_Sva))/COD_feed;
            j=j+1;
        end
        i=i+1;
    end
    
     for i=2:length(substrates_ratio)-1
        for j=2:length(pH)-1
            VFAyield_matrix(i,j)=VFAyield_matrix(i,j)*k+(VFAyield_matrix(i+1,j)+VFAyield_matrix(i-1,j)+VFAyield_matrix(i,j+1)+VFAyield_matrix(i,j-1))*(1-k)/4;
            j=j+1;
        end
        i=i+1;
     end
    
     figure(fig_number);
    f_my_contorno2(VFAyield_matrix,substrates_ratio,pH)
    title('(b)')
    i=i+1;
else
end


%% Acidification map

if flag3==1
    i=1;
    for i=1:length(substrates_ratio)
        j=1;
        for j=1:length(pH)
            aux1=sum(states_matrix(j+length(pH)*(i-1),pos_Ssu:pos_Sfa));
            aux2=sum(states_matrix(j+length(pH)*(i-1),pos_Xc1:pos_Xc2));
            aux3=sum(states_matrix(j+length(pH)*(i-1),pos_Xch:pos_Xi));
            acidification_matrix(j,i)=100*(COD_feed-aux1-aux2-aux3)/COD_feed;
            j=j+1;
        end
        i=i+1;
    end
    
     for i=2:length(substrates_ratio)-1
        for j=2:length(pH)-1
            acidification_matrix(i,j)=acidification_matrix(i,j)*k+(acidification_matrix(i+1,j)+acidification_matrix(i-1,j)+acidification_matrix(i,j+1)+acidification_matrix(i,j-1))*(1-k)/4;
            j=j+1;
        end
        i=i+1;
    end
    
      figure (fig_number+1)
    f_my_contorno2(acidification_matrix,substrates_ratio,pH)
    title('(a)')
    i=i+1;
else
end


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
    for i=1:length(substrates_ratio)
        j=1;
        for j=1:length(pH)
            Sac=states_matrix(j+length(pH)*(i-1),pos_Sac)/64;
            Spro=states_matrix(j+length(pH)*(i-1),pos_Sac+1)/112;
            Sbu=states_matrix(j+length(pH)*(i-1),pos_Sac+2)/160;
            Sva=states_matrix(j+length(pH)*(i-1),pos_Sva)/208;
            if Sac>Spro
                HV=Sva+Spro;
                HB=Sbu+(Sac-Spro)/2;
            else
                HV=Sva+Sac;
                HB=Sbu;
            end
            spectrum_matrix(j,i)=HV/(HB+HV);
            j=j+1;
        end
        i=i+1;
    end
    
    for i=2:length(substrates_ratio)-1
        for j=2:length(pH)-1
            spectrum_matrix(i,j)=spectrum_matrix(i,j)*k+(spectrum_matrix(i+1,j)+spectrum_matrix(i-1,j)+spectrum_matrix(i,j+1)+spectrum_matrix(i,j-1))*(1-k)/4;
            j=j+1;
        end
        i=i+1;
    end
    
    handle =figure(fig_number+3)
    f_my_contorno2(spectrum_matrix,substrates_ratio,pH)
    title('(d)')
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
    for i=1:length(substrates_ratio)
        j=1;
        for j=1:length(pH)
            Qliq=Qliq_vector(j+length(pH)*(i-1));
            productivity_matrix(j,i)=sum(states_matrix(j+length(pH)*(i-1),pos_Sac:pos_Sva))*Qliq/Vliq;
            j=j+1;
        end
        i=i+1;
    end
    
    for i=2:length(substrates_ratio)-1
        for j=2:length(pH)-1
            productivity_matrix(i,j)=productivity_matrix(i,j)*k+(productivity_matrix(i+1,j)+productivity_matrix(i-1,j)+productivity_matrix(i,j+1)+productivity_matrix(i,j-1))*(1-k)/4;
            j=j+1;
        end
        i=i+1;
    end
    
    figure(fig_number+4)
    f_my_contorno2(productivity_matrix,substrates_ratio,pH)
    title('(c)')
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
        for j=1:length(substrates_ratio)
            k=1;
            for k=1:length(pH)
                flag5_matrix(k,j)=inhibition_matrix(k+length(pH)*(j-1),flag5_pos);
                k=k+1;
            end
            j=j+1;
        end
        fig_number=i+6;
        figure (fig_number)
        f_my_contorno2(flag5_matrix,substrates_ratio,pH)
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
        for j=1:length(substrates_ratio)
            k=1;
            for k=1:length(pH)
                flag1_matrix(k,j)=states_matrix(k+length(pH)*(j-1),flag1_pos);
                k=k+1;
            end
            j=j+1;
        end
        fig_number=length(flag5)+i+6;
       figure (fig_number)
        f_my_contorno2(flag1_matrix,substrates_ratio,pH)
        figure_name=Abb(flag1_pos);
        title(figure_name);
        i=i+1;
    end
end


end
