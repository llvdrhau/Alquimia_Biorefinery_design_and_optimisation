% loadModelXls.m -Load all the parameters from Excel sheet
function [R] = loadDataLCAXls(idR)
warning off all

% If only one reactor is used is set to '' and no extra for the names of the xls sheets is needed
if idR==1,  id = '';    else    id = num2str(idR);   end
route = char('ResultadosLCA.xlsx');

% GENERAL MODEL PARAMETERS from the Excel file
fprintf('\n> Loading and creating GENERAL MODEL PARAMETERS...')

[LCAResults, pLCARNames]   = xlsread(route)
aux = '';