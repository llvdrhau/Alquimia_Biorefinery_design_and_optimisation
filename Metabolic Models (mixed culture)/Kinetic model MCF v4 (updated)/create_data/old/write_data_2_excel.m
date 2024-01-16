%% training data 
load('Data_4_ML_model')

% input data
T = array2table(inputMatrix);
T.Properties.VariableNames(1:size(inputMatrix,2)) = nameInputs;

% output data
T2 = array2table(outputMatrix);
T2.Properties.VariableNames(1:size(outputMatrix,2)) = nameOutputs;

% export to excel 
writetable(T,'TrainingData.xlsx','Sheet','input');
writetable(T2,'TrainingData.xlsx','Sheet','output');

%% validation data 

load('Validation_Data_4_ML_model.mat')
T = array2table(inputMatrix);
T.Properties.VariableNames(1:size(inputMatrix,2)) = nameInputs;

% output data
T2 = array2table(outputMatrix);
T2.Properties.VariableNames(1:size(outputMatrix,2)) = nameOutputs;

% export to excel 
writetable(T,'ValidationData.xlsx','Sheet','input');
writetable(T2,'ValidationData.xlsx','Sheet','output');