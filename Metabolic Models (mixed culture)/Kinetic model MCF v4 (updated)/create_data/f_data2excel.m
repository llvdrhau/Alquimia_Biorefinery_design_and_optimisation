function [] = f_data2excel(data, colNames, excelName, sheetName)
%UNTITLED Summary of this function goes here
%   turns data into excel file 

T = array2table(data);
T.Properties.VariableNames(1:size(data,2)) = colNames;
writetable(T,excelName,'Sheet',sheetName);

end