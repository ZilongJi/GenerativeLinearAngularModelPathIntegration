%% FittingData 
% This script will lo%% FittingData 
% This script will load the data, rotate the paths and call the solver to find the parameters 
% that minimize the mean distance to generated x3' which matchs the physical leg 3

%% Cleaning variables
clearvars; clear all; close all; clc;

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
%load('Data/AllDataErrors2018_V2.mat');
%V3 including the missing participants, which should always being loaded now
load('Data/AllDataErrors2018_V3.mat');
savefolder = "C:/Users/Zilong/Desktop/path integration model/Andrea's matlab code/GammaModelAllReplaceWithBias/Output/";

%% setting the configuration
config.UseGlobalSearch = true;
config.USEOoBTrials = true;
%load the data
if config.USEOoBTrials== true 
    YoungControls = TransformPathsOoB(YoungControls); 
else
    YoungControls = TransformPaths(YoungControls);
end

resultfolder = savefolder+"ModelComp/";
config.ResultFolder = resultfolder;
%create storing folder for trajectory if not exist
if ~exist(resultfolder, 'dir')
   mkdir(resultfolder);
end

%% the allocentric PI model without weber's law
config.ModelName = "allo";
config.NumParams = 1;
[~, ~, ~, ~, AllYoungIC_Allocentric] = getResultsAllConditions(YoungControls, config);

%% the allocentric PI model with weber's law
config.ModelName = "allo_weber";
config.NumParams = 1;
[~, ~, ~, ~, AllYoungIC_AllocentricWeber] = getResultsAllConditions(YoungControls, config);

%% the egocentric PI model
config.ModelName = "ego";
config.NumParams = 2;

%% Our base model
%no consider of the rotation gain factor in leg2, i.e., gamma, G3=1, g2=1, g3, k3, sigma, nu. #params=6
config.ModelName = "set_g2_1";
config.NumParams = 5;
[~, ~, ~, ~, AllYoungIC_Base] = getResultsAllConditions(YoungControls, config);

%% Box Plot of AIC
ICType = "AIC";

IC_Allocentric = getRawIC(AllYoungIC_Allocentric, ICType);
IC_AllocentricWeber = getRawIC(AllYoungIC_AllocentricWeber, ICType);
IC_Base = getRawIC(AllYoungIC_Base, ICType);

IC_Cond1 = [IC_Allocentric{1}',IC_AllocentricWeber{1}',IC_Base{1}'];
plotBoxPlot(IC_Cond1, "NoChange", ICType, resultfolder);

IC_Cond2 = [IC_Allocentric{2}',IC_AllocentricWeber{2}',IC_Base{2}'];
plotBoxPlot(IC_Cond2, "NoDistalCue", ICType, resultfolder);

IC_Cond3 = [IC_Allocentric{3}',IC_AllocentricWeber{3}',IC_Base{3}'];
plotBoxPlot(IC_Cond3, "NoOpticalFlow", ICType, resultfolder);

% Box Plot of BIC 
ICType = "BIC";
label_xaxis_data = ['allo', 'allo_weber', 'base'];

IC_Allocentric = getRawIC(AllYoungIC_Allocentric, ICType);
IC_AllocentricWeber = getRawIC(AllYoungIC_AllocentricWeber, ICType);
IC_Base = getRawIC(AllYoungIC_Base, ICType);

IC_Cond1 = [IC_Allocentric{1}',IC_AllocentricWeber{1}',IC_Base{1}'];
plotBoxPlot(IC_Cond1, "NoChange", ICType, resultfolder);

IC_Cond2 = [IC_Allocentric{2}',IC_AllocentricWeber{2}',IC_Base{2}'];
plotBoxPlot(IC_Cond2, "NoDistalCue", ICType, resultfolder);

IC_Cond3 = [IC_Allocentric{3}',IC_AllocentricWeber{3}',IC_Base{3}'];
plotBoxPlot(IC_Cond3, "NoOpticalFlow", ICType, resultfolder);

% Box Plot of NegLogLikelihood
ICType = "NegLogLikelihood";

IC_Allocentric = getRawIC(AllYoungIC_Allocentric, ICType);
IC_AllocentricWeber = getRawIC(AllYoungIC_AllocentricWeber, ICType);
IC_Base = getRawIC(AllYoungIC_Base, ICType);

IC_Cond1 = [IC_Allocentric{1}',IC_AllocentricWeber{1}',IC_Base{1}'];
plotBoxPlot(IC_Cond1, "NoChange", ICType, resultfolder);

IC_Cond2 = [IC_Allocentric{2}',IC_AllocentricWeber{2}',IC_Base{2}'];
plotBoxPlot(IC_Cond2, "NoDistalCue", ICType, resultfolder);

IC_Cond3 = [IC_Allocentric{3}',IC_AllocentricWeber{3}',IC_Base{3}'];
plotBoxPlot(IC_Cond3, "NoOpticalFlow", ICType, resultfolder);

%% A function for getting Results from All Conditions
function [AllParams, AllX, AllDX, AllTheta, AllIC] = getResultsAllConditions(TransformedData, config)
    %get the estimated parameters, X, DX, Theta, IC for all trial
    %conditions for each group of data.
    AllParams = cell(0); AllX = cell(0);AllDX = cell(0); AllTheta = cell(0); AllIC = cell(0);
    for TRIAL_FILTER=1:3
        config.TrialFilter = TRIAL_FILTER;
        tic
        disp('%%%%%%%%%%%%%%% PERFORMING FITTING %%%%%%%%%%%%%%%');
        Results = PerformGroupFit(TransformedData, config);
        AllParams{TRIAL_FILTER} = Results.estimatedParams; 
        AllX{TRIAL_FILTER} = Results.X;
        AllDX{TRIAL_FILTER}=Results.DX;      
        AllTheta{TRIAL_FILTER}=Results.THETADX;
        AllIC{TRIAL_FILTER}=Results.IC;
        toc
    end
end    

%% function for get the IC data out from the data structure
function IC_AllConds=getRawIC(IC, ICType)
    IC_AllConds = cell(1,3);
    for i=1:3 %three conditions
        if ICType == "AIC"
            icdata = [];
            for j=1:length(IC{i})
                icdata = [icdata,IC{i}{j}.aic];
            end
        elseif ICType == "BIC"
            icdata = [];
            for j=1:length(IC{i})
                icdata = [icdata,IC{i}{j}.bic];
            end
        elseif ICType == "NegLogLikelihood"
            icdata = [];
            for j=1:length(IC{i})
                icdata = [icdata,IC{i}{j}.negll];
            end      
        else
            error("Choose correct IC type!")
        end
        IC_AllConds{i}=icdata;
    end
end

%% function for Box plot
function plotBoxPlot(data, CondType, ICType, resultfolder)

    f = figure('visible','off','Position', [100 100 500 300]);
    
    %%%set paramsters
    num_boxplots = size(data,2);
    box_lineWidth = 0.5;
    whisker_value = 1.5;
    box_widths_value = 0.3;
    box_color_transparency = 0.5; %faceAlpha
    median_lineWidth = 2;
    median_color = 'k';
    scatter_jitter_value = 0.1;
    scatter_markerSize=10;
    scatter_marker_edgeColor = 'k';
    scatter_marker_edgeWidth = 0.5;
    scatter_color_transparency = 0.7; %faceAlpha
%     outlier_marker = 's';
%     outlier_jitter_value = 0.1;
%     outlier_marker_edgeColor = 'k';
%     outlier_markerSize = 4; 
%     outlier_marker_edgeWidth = 0.5;

    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)

    %%% Color definition %%%
    color_scheme_npg = [0.9020    0.2941    0.2078; ...
                        0.3020    0.7333    0.8353; ...
                        0         0.6275    0.5294; ...
                        0.2353    0.3294    0.5333; ...
                        0.9529    0.6078    0.4980; ...
                        0.5176    0.5686    0.7059; ...
                        0.5686    0.8196    0.7608; ...
                        0.8627         0         0; ...
                        0.4941    0.3804    0.2824; ...
                        0.6902    0.6118    0.5216 ];    
    box_colors = color_scheme_npg(1:num_boxplots,:);
    
    %% main boxplot one box for each column in data
    bp = boxplot(data, 'whisker',whisker_value,'symbol','', ... %symbol ='' making outlier invisible
        'color','k', 'labels', {'allo', 'allo_weber', 'base'},...
        'widths',box_widths_value);
    
    set(bp,'linewidth',box_lineWidth);

    if ICType=="AIC"
        ylabel('Akaike information criterion (AIC)');
    elseif ICType=="BIC"
        ylabel('Bayesian Inference Criterion (BIC)');
    elseif ICType=="NegLogLikelihood"
        ylabel('Negative Loglikelihood');
    else
        error("Choose correct IC type!");
    end
    
    %% Coloring each box
    h = findobj(gca,'Tag','Box');
    for i = 1:length(h)
        %note that first getting the most right box, so do "length(h)-i+1"
        patch(get(h(i),'XData'),get(h(i),'YData'),box_colors(length(h)-i+1,:),'FaceAlpha',box_color_transparency);
    end
    % Sending patch to back of the figure so that median can be drawn on top of it
    set(gca,'children',flipud(get(gca,'children'))) 

    %% Adjusting median
    h=findobj(gca,'tag','Median');
    for i = 1:length(h)
        h(i).LineWidth = median_lineWidth;
        h(i).Color = median_color;
    end

    %% add scatter plot of the data
    num_points = size(data,1);
    for i=1:size(data,2)
        hold on
        x = i*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, data(:,i), scatter_markerSize, ...
                'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',box_colors(i,:), ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 
    end

%     %% Adjusting outliers
%     h=findobj(gca,'tag','Outliers');
%     for i = 1:length(h)
%         h(i).MarkerFaceColor = [0.5,0.5,0.5]; %alpha(0.7)
%         h(i).MarkerEdgeColor = outlier_marker_edgeColor;
%         h(i).MarkerSize = outlier_markerSize;
%         h(i).LineWidth = outlier_marker_edgeWidth;
%     end

    %% Further post-processing the figure
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , 1:1:100,... 
        'LineWidth'   , .5        );

    %% save figure
    exportgraphics(f,resultfolder+"/boxplot_"+ICType+CondType+".png",'Resolution',300);
end
