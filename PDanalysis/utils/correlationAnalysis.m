function [V,VPALL,Tcomf]=correlationAnalysis(T,dataPosition,methodsel,scales,scaleType)
% correlationAnalysis - Computes the correlation between activity features and clinical scales.
%
% Description:
% This function performs correlation analysis between extracted activity features
% (e.g., derived from trunk or wrist sensors) and clinical scales for Parkinson’s Disease (PD).
% It supports various statistical methods such as Pearson, Spearman, Kendall, regression, and linear mixed-effects models.
% The results include correlation coefficients and p-values, which are saved to files for further analysis.
%
% Inputs:
% - T: (table) A table containing activity features. The first column is `DataID`, and subsequent
%      columns contain extracted activity features for each participant.
% - dataPosition: (string) Indicates the sensor position ('Body' or 'Wrist') for analysis.
% - methodsel: (string) The statistical method to use for correlation analysis. Supported methods:
%   - "Pearson": Pearson correlation.
%   - "Spearman": Spearman rank correlation.
%   - "Kendall": Kendall rank correlation.
%   - "regress": Linear regression.
%   - "fitlme": Linear mixed-effects model.
% - scales: (struct) A structure containing clinical scales for each participant. The field names
%           should match the participant IDs, and each field should contain a vector of scale values.
%
% Outputs:
% - V: (table) A table of correlation coefficients for each activity feature and clinical scale.
% - VPALL: (table) A table combining correlation coefficients and p-values in the format `value(pvalue)`.
% - Tcomf: (table) A table containing extracted features for each participant averaged across 24h chunks

% Key Features:
% 1. Computes correlation coefficients and p-values for each activity feature and clinical scale.
% 2. Supports advanced methods like linear mixed-effects models (`fitlme`) for hierarchical data.
% 3. Formats results as both raw numerical values (`V`) and human-readable strings (`P`).
% 4. Automatically saves results to `.mat` and `.xlsx` files in the `results/` directory.

% Workflow:
% 1. Extracts activity features (`fcorrALL`) for each participant, aggregating multiple samples by averaging.
% 2. Matches activity features with clinical scales (`scales`).
% 3. Computes correlation coefficients and p-values for all activity feature-clinical scale pairs.
% 4. Saves the results as:
%    - A `.mat` file containing `values` (correlation coefficients) and `pvalues`.
%    - An `.xlsx` file with human-readable formatted results (e.g., `value(pvalue)`).

% Example Usage:
% [V, P] = correlationAnalysis(T, 'Body', 'Pearson', scales);
% This computes Pearson correlations between activity features from the body sensor
% and clinical scales stored in the `scales` structure.

DataID=table2array(T(:,1));
fcorrALL=table2array(T(:,3:end));
ParID=unique(DataID);
fcorrALL1=zeros(length(ParID),size(fcorrALL,2));

scalesNames=scales.Scale;

for i=1:length(ParID)
    % for i=2:2
    sel=DataID==ParID(i);
    if sum(sel)==1
        fcorrALL1(i,:)=fcorrALL(sel,:);

    else
        if dataPosition=="Body"
            fcorrALL1(i,:)=mean(fcorrALL(sel,:));


        elseif dataPosition=="Wrist"
            fcorrALL1(i,:)=mean(fcorrALL(sel,:),"omitnan");

        end
    end
    disp(strcat(strcat(ParID(i),':',num2str(sum(sel)))))
    if i==1
        SL=scales.(ParID(i));
    else
        SL=[SL,scales.(ParID(i))];
    end

end

if dataPosition=="Body"
    load("Info\ActivityFeatureName.mat",'BodyfeatureName')
    featureName=BodyfeatureName;
elseif dataPosition=="Wrist"
    load("Info\ActivityFeatureName.mat",'WristfeatureName')
    featureName=WristfeatureName;
end
% Create the main table for fcorrALL
Tcomf = array2table(fcorrALL1, 'VariableNames', cellstr(featureName));

% Add DataID as the first column
Tcomf = [table(ParID, 'VariableNames', {'DataID'}), Tcomf];

SL=SL';

for j=1:size(SL,2)
    a=SL(:,j);
    if methodsel=="fitlme"
        stat=[0,0,0,0];
    else
        stat=[0,0];
    end

    for i=1:size(fcorrALL1,2)
        b=fcorrALL1(:,i)  ;
        if methodsel=="Pearson"
            % Pearson corre
            [R,P] = corrcoef(a,b);
            stat=[stat;R(1,2),P(1,2)];
        elseif or(methodsel=="Spearman",methodsel=="Kendall")
            [r,p]=corr(a,b,"type",methodsel);
            stat=[stat;r,p];
        elseif methodsel=="regress"
            if ismember(i,65:70)
                stat=[stat;[nan,nan]];
            else
                [~,~,~,~,stats] = regress(b,[a,ones(size(a))]);
                stat=[stat;[stats(1),stats(3)]];
            end
        elseif methodsel=="fitlme"
            data = table(DataID(rangesel), b, a);
            % stat=[stat;[ model.Coefficients.Estimate(2), model.Coefficients.pValue(2),stats{1}.Estimate]];
            lme = fitlme(data, 'b ~ a + (1|Var1) ');
            [~,~,stats]= fixedEffects(lme);
            randomEffectsTable = randomEffects(lme);
            % Extract coefficients and p-values
            stat=[stat;[ lme.Coefficients.Estimate(2), lme.Coefficients.SE(2),lme.Coefficients.tStat(2),lme.Coefficients.pValue(2)]];
        end

    end
    stat(1,:)=[];
    if j==1
        STAT=stat;
    else
        STAT=[STAT,stat];
    end

end
values=STAT(:,1:2:end);
pvalues=STAT(:,2:2:end);
if dataPosition=="Body"
    load("Info\ActivityFeatureName.mat",'BodyfeatureName')
    featureName=BodyfeatureName;
elseif dataPosition=="Wrist"
    load("Info\ActivityFeatureName.mat",'WristfeatureName')
    featureName=WristfeatureName;
end

rej=any(isnan(values));
scalesNames=scalesNames(~rej);
values=values(:,~rej);
pvalues=pvalues(:,~rej);

% Combine values and pvalues into formatted strings
formattedStrings = arrayfun(@(v, p) sprintf('%.2f(%.3f)', v, p), values, pvalues, 'UniformOutput', false);

% Create the main table for fcorrALL
V = array2table(values, 'VariableNames',  cellstr(scalesNames));
V = [table(featureName', 'VariableNames', {'ActivityFeature'}),  V];

VPALL = array2table(formattedStrings, 'VariableNames',  cellstr(scalesNames));
VPALL = [table(featureName', 'VariableNames', {'ActivityFeature'}), VPALL];

savepath=fullfile("results",strcat(dataPosition,methodsel,'_',scaleType));
save(savepath,"pvalues","values")
fprintf('Saved correlation results in %s.\n', strcat(savepath,'.mat'));
% Save as CSV (converting cell array to table if needed)
% Write the table to an Excel file
writetable(VPALL , strcat(savepath,'.xlsx'));
fprintf('Saved correlation results in %s.\n', strcat(savepath,'.xlsx'));
figure_name=strcat(savepath,'.png');
correlationPlot(values,pvalues,figure_name,'labels',scalesNames,featureName)

end