function stoiMatrix = f_my_stoichiometric_matrix(stoiTable, stoiPar, stoiNames)

% This function builds the stoichiometric matrix from the stoichiometric
% parameters and the stoichiometric table
% *********************************************************************** %
% INPUTS:
% stoiTable 
% stoiPar
% *********************************************************************** %
% OUTPUT:
% stoiMatrix
% *********************************************************************** %

% Mateo Saavedra del Oso. CRETUS Institute. University of Santiago de Compostela.
% Spain. November 2020. Please contact msaavedra.deloso@usc.es if you intend 
% to use this code. 


%% loading stoichiometric parameters

Fch_xc1 = stoiPar(strcmp(stoiNames,'Fch_xc1'));     % carbohydrates fraction in composite 1
Fpr_xc1 = stoiPar(strcmp(stoiNames,'Fpr_xc1'));     % proteins fraction in composite 1
Ftdf_xc1 = stoiPar(strcmp(stoiNames,'Ftdf_xc1'));   % TDF fraction in composite 1
Fli_xc1 = stoiPar(strcmp(stoiNames,'Fli_xc1'));     % lipids fraction in composite 1
Fxi_xc1 = stoiPar(strcmp(stoiNames,'Fxi_xc1'));     % inerts fraction in composite 1
Fsi_xc1 = stoiPar(strcmp(stoiNames,'Fsi_xc1'));     % soluble inerts fraction in composite 1

Fch_xc2 = stoiPar(strcmp(stoiNames,'Fch_xc2'));     % carbohydrates fraction in composite 2
Fpr_xc2 = stoiPar(strcmp(stoiNames,'Fpr_xc2'));     % proteins fraction in composite 2
Ftdf_xc2 = stoiPar(strcmp(stoiNames,'Ftdf_xc2'));   % TDF fraction in composite 2
Fli_xc2 = stoiPar(strcmp(stoiNames,'Fli_xc2'));     % lipids fraction in composite 2
Fxi_xc2 = stoiPar(strcmp(stoiNames,'Fxi_xc2'));     % inerts fraction in composite 2
Fsi_xc2 = stoiPar(strcmp(stoiNames,'Fsi_xc2'));     % soluble inerts fraction in composite 2

Fch_xc3 = stoiPar(strcmp(stoiNames,'Fch_xc3'));     % carbohydrates fraction in composite 3
Fpr_xc3 = stoiPar(strcmp(stoiNames,'Fpr_xc3'));     % proteins fraction in composite 3
Ftdf_xc3 = stoiPar(strcmp(stoiNames,'Ftdf_xc3'));   % TDF fraction in composite 3
Fli_xc3 = stoiPar(strcmp(stoiNames,'Fli_xc3'));     % lipids fraction in composite 3
Fxi_xc3 = stoiPar(strcmp(stoiNames,'Fxi_xc3'));     % inerts fraction in composite 3
Fsi_xc3 = stoiPar(strcmp(stoiNames,'Fsi_xc3'));     % soluble inerts fraction in composite 3

Ffa_li = stoiPar(strcmp(stoiNames,'Ffa_li'));       % lipids fraction converted into LCFA
Fet_af = stoiPar(strcmp(stoiNames,'Fet_af'));       % ethanol fraction from acidogenic fermentation
Fva_af = stoiPar(strcmp(stoiNames,'Fva_af'));       % valerate fraction from acidogenic fermentation
Fbu_af = stoiPar(strcmp(stoiNames,'Fbu_af'));       % butyrate fraction from acidogenic fermentation
Fpro_af = stoiPar(strcmp(stoiNames,'Fpro_af'));     % propionate fraction from acidogenic fermentation
Fac_af = stoiPar(strcmp(stoiNames,'Fac_af'));       % acetate fraction from acidogenic fermentation
Fh2_af = stoiPar(strcmp(stoiNames,'Fh2_af'));       % hydrogen fraction from acidogenic fermentation

Yaf = stoiPar(strcmp(stoiNames,'Yaf'));             % yield acidogenic fermentation biomass
Yfa = stoiPar(strcmp(stoiNames,'Yfa'));             % yield LCFA degraders
Yc4 = stoiPar(strcmp(stoiNames,'Yc4'));             % yield butyrate and valerate degraders
Ypro = stoiPar(strcmp(stoiNames,'Ypro'));           % yield propionate degraders
Yac = stoiPar(strcmp(stoiNames,'Yac'));             % yield acetate degraders
Yh2 = stoiPar(strcmp(stoiNames,'Yh2'));             % yield hydrogen degraders

Cxc1 = stoiPar(strcmp(stoiNames,'Cxc1'));                  % total inorganic carbon in composite 1
Cxc2 = stoiPar(strcmp(stoiNames,'Cxc2'));                  % total inorganic carbon in composite 2
Cxc3 = stoiPar(strcmp(stoiNames,'Cxc3'));                  % total inorganic carbon in composite 3
Cbac = stoiPar(strcmp(stoiNames,'Cbac'));                  % total inorganic carbon in active biomass
Csi = stoiPar(strcmp(stoiNames,'Csi'));                    % total inorganic carbon in soluble inorganics
Cxi = stoiPar(strcmp(stoiNames,'Cxi'));                    % total inorganic carbon in particulate inorganics
Cch = stoiPar(strcmp(stoiNames,'Cch'));                    % total inorganic carbon in carbohydrates
Cpr=stoiPar(strcmp(stoiNames,'Cpr'));                    % total inorganic carbon in proteins
Cli=stoiPar(strcmp(stoiNames,'Cli'));                    % total inorganic carbon in lipids
Csu=stoiPar(strcmp(stoiNames,'Csu'));                    % total inorganic carbon in sugars
Caa=stoiPar(strcmp(stoiNames,'Caa'));                    % total inorganic carbon in aminoacids
Cfa=stoiPar(strcmp(stoiNames,'Cfa'));                    % total inorganic carbon in LCFA
Cva=stoiPar(strcmp(stoiNames,'Cva'));                    % total inorganic carbon in valerate
Cpro=stoiPar(strcmp(stoiNames,'Cpro'));                  % total inorganic carbon in propionate
Cbu=stoiPar(strcmp(stoiNames,'Cbu'));                    % total inorganic carbon in butyrate
Cac=stoiPar(strcmp(stoiNames,'Cac'));                    % total inorganic carbon in acetate
Cch4=stoiPar(strcmp(stoiNames,'Cch4'));                  % total inorganic carbon in methane

Nxc1=stoiPar(strcmp(stoiNames,'Nxc1'));                  % total inorganic nitrogen in composite 1
Nxc2=stoiPar(strcmp(stoiNames,'Nxc2'));                  % total inorganic nitrogen in composite 2
Nxc3=stoiPar(strcmp(stoiNames,'Nxc3'));                  % total inorganic nitrogen in composite 3
Nbac=stoiPar(strcmp(stoiNames,'Nbac'));                  % total inorganic nitrogen in active biomass
NI=stoiPar(strcmp(stoiNames,'NI'));                      % total inorganic nitrogen in inorganics
Naa=stoiPar(strcmp(stoiNames,'Naa'));                    % total inorganic nitrogen in aminoacids

%% building stoichiometric matrix

stoiTable=stoiTable(3:(end-1),3:(end-1));
[np,nc]=size(stoiTable);
stoiMatrix=zeros(np,nc);

for i=1:np
    for j=1:nc
        elementCell=stoiTable{i,j};
        if isa(elementCell,'double')
            stoiMatrix(i,j) = elementCell;
        else
            stoiMatrix(i,j)=eval(elementCell);
        end
    end
end

end