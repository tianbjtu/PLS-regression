% clear;clc;
%%
DATA_DIR ='E:\code\data\';
FOLD_NUM = 10; DIM = 50;
load (strcat (DATA_DIR, 'HCP_pCorr_1003s_200ROI'));
load (strcat (DATA_DIR, 'HCP_1003s_fMRI_ID'));

%%
%-------------------Age+Edu+3IQ
load (strcat (DATA_DIR, 'HCP_SubjInfo_Tian_1206_45Items'));AA = SubjInfo(:,[1, 3, 6]); 
load (strcat (DATA_DIR, 'HCP_Cognition_Tian_1206_83Items')); BB=SubjInfo(:,[46, 50, 52]);
SubjInfo=[AA,BB]; 
%% 
% Delete samples whose fMRI data unavailable 
ItemNum = size(SubjInfo, 2);
for Tmp = 1 : size(SubjInfo, 1)
    if(find(fMRI_ID == SubjInfo(Tmp, 1)))
        SubjInfo(Tmp, ItemNum+1) = 1;
    else
        SubjInfo(Tmp, ItemNum+1) = 0;
    end
end
SubjNoUseTmp = find((SubjInfo(:, end) == 1)');   SubjInfoP = SubjInfo(SubjNoUseTmp, 1:end-1); %fMRI data available
% Delete Saples whose Labels unavailable
SubjNoUse = find (sum(SubjInfoP'==-9999) == 0); 
SampleSize = length(SubjNoUse);


%%
CorrArray = CorrVec(SubjNoUse, :);
% Label = [SubjInfoP(SubjNoUse, [2:end])]; %% 5_label
Label_zong = [SubjInfoP(SubjNoUse, [2:end])];
% Label = [SubjInfoP(SubjNoUse, [2])]; %% age label
% Label = [SubjInfoP(SubjNoUse, [3])]; %% Edu
% Label = [SubjInfoP(SubjNoUse, [4])]; %% Fluid
% Label = [SubjInfoP(SubjNoUse, [5])]; %% all_IQ
Label = [SubjInfoP(SubjNoUse, [6])]; %% Crystal
%%
[TrainNo, TestNo, Permed] = f_NFold(SampleSize, FOLD_NUM);
Predicted = []; LabelPermed = [];
bz=[];
bsingle = [];
for FoldNo = 1 : FOLD_NUM
    FoldNo
    LabelTrain = Label(TrainNo{FoldNo}, :);
    CorrTrain = CorrArray(TrainNo{FoldNo}, :);
    Train_Num = length(TrainNo{FoldNo});
    Test_Num = length(TestNo{FoldNo});
    
    LabelMean = mean(LabelTrain, 1); 
    CorrMean = mean(CorrTrain, 1);

    [XLOADINGS,YLOADINGS,XSCORES,YSCORES, Beta,pctVar] = plsregress(CorrTrain, LabelTrain, DIM);
    
    CorrTest = CorrArray(TestNo{FoldNo}, :); 
    CorrTest_Norm = CorrTest - repmat(CorrMean, size(CorrTest, 1), 1);    
    B = Beta(2:end, :);
    TestLabel_Pred = CorrTest_Norm * B;
    
    %%%% ===== b
    bz = [bz,B];
    
    Predicted = [Predicted; TestLabel_Pred + repmat(LabelMean, size(TestLabel_Pred, 1), 1)];
    LabelPermed = [LabelPermed; Label(TestNo{FoldNo}, :);];
end
 r = corrcoef(Predicted, LabelPermed)


%%
for Tmp = 1 : size(LabelPermed, 2)
    %%%----------Corr£¬RMSE£¬CoD
    r = corrcoef(Predicted(:, Tmp), LabelPermed(:, Tmp));
    R(Tmp) = r(1,2); 
    RMSE(Tmp) = sqrt (sum ((Predicted(:, Tmp) - LabelPermed(:, Tmp)) .^ 2) / (length(Label)-1));
    mm = mean(LabelPermed(:, Tmp));
    sse = sum((LabelPermed(:, Tmp)-Predicted(:, Tmp)).^2);
    sst = sum((LabelPermed(:, Tmp)-mm).^2); 
    CoD(Tmp) = 1-sse/sst;
end

R(1:end)
CoD(1:end)
