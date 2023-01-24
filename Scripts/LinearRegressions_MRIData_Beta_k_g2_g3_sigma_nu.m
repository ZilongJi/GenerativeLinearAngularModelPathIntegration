%% Preparing the data
VAM_PrepareBaseConfig

%% Preprocessing the data
VAM_PreprocessData

%% Setting the model we are interested in
% Eventually modify config paramteters we are interested in. For example
% for this graph we are not interested in running the model with splitted
% conditions so we will set the relative config to false 
% force it to not run
rng("default")
config.useTrialFilter = true;
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName); %length(config.ParamName); % Set 100 here to avoid producing the model
% Run the model
VAM

%% Model run completed, preparing the data for plotting figures
config.ResultFolder = pwd + "/Output/MRIvsModelParamsLinearRegressions/Beta_k_g2_g3_sigma_nu";
% Create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

ColorPattern;

%% Getting Information from results:
YoungControlsParameters   = averageAcrossConditions(YoungControls.Results.estimatedParams);
HealthyControlsParameters = averageAcrossConditions(HealthyControls.Results.estimatedParams);
MCIUnkParameters          = averageAcrossConditions(MCIUnk.Results.estimatedParams);
MCINegParameters          = averageAcrossConditions(MCINeg.Results.estimatedParams);
MCIPosParameters          = averageAcrossConditions(MCIPos.Results.estimatedParams);

%% Preparing linear mixed effect table
HealthyControls_lmetable = createLMETable(HealthyControls.MRI,HealthyControlsParameters);
MCIUnk_lmetable          = createLMETable(MCIUnk.MRI,MCIUnkParameters);
MCINeg_lmetable          = createLMETable(MCINeg.MRI,MCINegParameters);
MCIPos_lmetable          = createLMETable(MCIPos.MRI,MCIPosParameters);

MRIModelParamsDataTable  = [HealthyControls_lmetable; MCIUnk_lmetable; MCINeg_lmetable; MCIPos_lmetable];

%Removing nans (either from Nan parameters model or missing mri data)
filter = isnan(MRIModelParamsDataTable.beta) | isnan(MRIModelParamsDataTable.MCI);
% Counting removed samples
disp("Total Removed samples = " + sum(filter));
MRIModelParamsDataTable = MRIModelParamsDataTable(~filter,:);

disp("Removed helthy controls = " + (height(HealthyControlsParameters) - sum(MRIModelParamsDataTable.CSF == "HC")) ); 
disp("Removed unknown = " + (height(MCIUnkParameters) - sum(MRIModelParamsDataTable.CSF == "Unknown")) ); 
disp("Removed positive = " + (height(MCIPosParameters) - sum(MRIModelParamsDataTable.CSF == "Positive")) ); 
disp("Removed negative = " + (height(MCINegParameters) - sum(MRIModelParamsDataTable.CSF == "Negative")) ); 

clear filter

%% Plotting linear regressions for ROIs
close all;
clc;

plotInfo.defaultTextSize = 14;
plotInfo.defaultLineSize = 2.0;
plotInfo.LineSizeMdl = 3.0;
plotInfo.titleFontSize = 20;
plotInfo.labelSize = 24;
plotInfo.axisSize = 25;
plotInfo.dataSize = 80;
plotInfo.legendFontSize = 18;
plotInfo.visible = "on";
plotInfo.ResultFolder = config.ResultFolder;
plotInfo.color_scheme_group = config.color_scheme_group;

% Getting the indices for the ROI
parameters_label = [{'\beta'}, {'k'}, {'g_2'}, {'g_3'}, {'\sigma'}, {'\nu'}];

mri_indeces = ["norm_ErC"...
    "norm_hippocampus" "norm_precuneus_volume_Dest"...
    "norm_inferiorparietal_volume_Dest" "norm_isthmuscingulate_volume_Dest"];
mri_indeces_label = ["Entorhinal Cortex"...
    "Hippocampus" "Precuneus"...
    "Inferior parietal" "Isthmus cingulate"];

addLine = [0, 1, 1, 1, 0, 0];

for param_i = 1:length(config.ParamName)    
    for mri_i = 1:length(mri_indeces)
        disp("Fitting " + parameters_label{param_i} + " vs " + mri_indeces_label(mri_i));

        pvalue = fitLinearRegression(...
            MRIModelParamsDataTable,...
            mri_indeces(mri_i),...
            mri_indeces_label(mri_i),...
            config.ParamName(param_i),...
            parameters_label{param_i},...
            addLine(param_i),...
            plotInfo...
        );
        plotLinearRegression(...
            MRIModelParamsDataTable,...
            mri_indeces(mri_i),...
            mri_indeces_label(mri_i),...
            config.ParamName(param_i),...
            parameters_label{param_i},...
            addLine(param_i),...
            plotInfo,...
            pvalue...
        );
    end
end
clear param_i mri_i

clear statuscode sex age

%% get the model parameters, average across the conditions
% remove nans if the row after mergin still contains nans
function dataout = averageAcrossConditions(data)

    dataout = [];
    pSize = length(data{1});
    paramsSize = width(data{1});

    for i = 1:pSize
        tempP = [];
        for j = 1:paramsSize
            tempP = [tempP mean([data{1}(i,j) data{2}(i,j) data{3}(i,j)],"omitnan")];
        end
        dataout = [dataout;tempP];
    end

end

%% Creates a unique table between mri and model parameters
function dataout = createLMETable(MRI, GroupParams)
    GroupParamsTable = array2table(GroupParams, "VariableNames", {'beta', 'k', 'g2', 'g3', 'sigma', 'nu'});
    dataout = [GroupParamsTable,MRI];
end

%%
function pvalue = fitLinearRegression(fullTable,x,xLabelPlot,y,yLabelPlot,addLine,plotInfo)

fullTable.Sex = nominal(fullTable.Sex);
fullTable.CSF = nominal(fullTable.CSF);
fullTable.MCI = nominal(fullTable.MCI);

fullTable.Age                               = zscore(fullTable.Age);  
fullTable.norm_ErC                          = zscore(fullTable.norm_ErC);
fullTable.norm_AntLat                       = zscore(fullTable.norm_AntLat);
fullTable.norm_PosMed                       = zscore(fullTable.norm_PosMed);
fullTable.norm_TE35                         = zscore(fullTable.norm_TE35);
fullTable.norm_hippocampus                  = zscore(fullTable.norm_hippocampus);
fullTable.norm_subiculum                    = zscore(fullTable.norm_subiculum);
fullTable.norm_posteriorCingulate           = zscore(fullTable.norm_posteriorCingulate);
fullTable.norm_inferiorparietal_volume_Dest = zscore(fullTable.norm_inferiorparietal_volume_Dest);
fullTable.norm_isthmuscingulate_volume_Dest = zscore(fullTable.norm_isthmuscingulate_volume_Dest);
fullTable.norm_superiorparietal_volume_Dest = zscore(fullTable.norm_superiorparietal_volume_Dest);
fullTable.norm_precuneus_volume_Dest        = zscore(fullTable.norm_precuneus_volume_Dest);

model_spec = y + " ~ Age + Sex + " + x;
mdl_full = fitlm(fullTable,model_spec) % To get the pvalue correcting for age and sex
anova_mdl_full = anova(mdl_full);
pvalue = anova_mdl_full.pValue(anova_mdl_full.Properties.RowNames == x);

if(pvalue < 0.05/5)
    disp("%%%%%%%%%%%%%%%%%%%%%%%");
    disp("%%%%%%%%%%%%%%%%%%%%%%% Survived correction " + xLabelPlot + " " + yLabelPlot + "%%%%%%%%%%%%%%%%%%%%%%%");
    disp("%%%%%%%%%%%%%%%%%%%%%%%");
end

end

%% Plotting a linear regression using a linear model
function plotLinearRegression(fullTable,x,xLabelPlot,y,yLabelPlot,addLine,plotInfo,pvalue)

    labels = ["HC", "Unknown", "Negative", "Positive"];
    fullTable.Sex = nominal(fullTable.Sex);
    xData = table2array(fullTable(:,find(strcmp(fullTable.Properties.VariableNames, x))));
    yData = table2array(fullTable(:,find(strcmp(fullTable.Properties.VariableNames, y))));
    csfStatus = fullTable.CSF;

    mdl = fitlm(xData,yData);

    % set figure info
    f = figure('visible', plotInfo.visible, 'Position', [0 0 500 400]);
    %%% Font type and size setting %%%
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',plotInfo.axisSize)
    set(0,'DefaultTextFontSize',plotInfo.axisSize)

    axSc = {};

    hold on;

    for i = 1:4
        axSc{i} = scatter(xData(csfStatus == labels(i)),yData((csfStatus == labels(i))), plotInfo.dataSize);
        axSc{i}.MarkerFaceColor = plotInfo.color_scheme_group(1+i,:);
        axSc{i}.MarkerEdgeColor = plotInfo.color_scheme_group(1+i,:) * 0.8;
        axSc{i}.MarkerFaceAlpha = 0.9; axSc{i}.MarkerEdgeAlpha = 0.9;
    end

    axLine = plot(mdl);
    axLine(1).Marker = 'none';
    axLine(2).Color = [0.2 0.2 0.2];
    axLine(3).Color = [0.2 0.2 0.2];
    axLine(4).Color = [0.2 0.2 0.2];
    axLine(2).LineWidth = plotInfo.LineSizeMdl;
    axLine(3).LineWidth = plotInfo.LineSizeMdl;
    axLine(4).LineWidth = plotInfo.LineSizeMdl;

    xlabel(xLabelPlot, Interpreter="tex");
    ylabel(yLabelPlot, Interpreter="tex");

    title("");
    % pvalue is Bonferroni corrected already
    if(pvalue < 0.001/5)
        title("***")
    elseif (pvalue < 0.01/5)
        title("**")
    elseif (pvalue < 0.05/5)
        title("*")
    end
    legplotInfo = legend([axSc{1}, axSc{2}, axSc{3}, axSc{4}], {'HC' 'MCI unk' 'MCI-' 'MCI+'}, "Location", "northeast", "AutoUpdate", "off");
    legplotInfo.FontSize = plotInfo.legendFontSize;
    
    ax = gca;

    if addLine > 0
        % Plotting y = 1
        tempX = [ax.XLim(1):0.0001:ax.XLim(2)];
        plot(tempX,addLine*ones(length(tempX)),"r--", LineWidth=plotInfo.LineSizeMdl);
    end

    ax.LineWidth = 3.0;
    ax.XAxis.FontSize = plotInfo.axisSize;
    ax.YAxis.FontSize = plotInfo.axisSize;
    ax.XLabel.FontSize = plotInfo.labelSize;
    ax.YLabel.FontSize = plotInfo.labelSize;
    ax.XAxis.ExponentMode = "manual";
    ax.XAxis.Exponent = -3;

    hold off;

    filename = convertCharsToStrings(yLabelPlot) + "vs" + convertCharsToStrings(xLabelPlot);
    
%     if exist(plotInfo.ResultFolder+"/"+filename+".png","file") == 2
%         delete(plotInfo.ResultFolder+"/"+filename+".png");
%     end

    exportgraphics(f,plotInfo.ResultFolder+"/"+filename+".png",'Resolution',300);
    exportgraphics(f,plotInfo.ResultFolder+"/"+filename+".pdf",'Resolution',300, 'ContentType','vector');

end