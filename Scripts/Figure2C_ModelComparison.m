%% Preparing the data
VAM_PrepareBaseConfig

%% Preprocessing the data
VAM_PreprocessData

%% 1, sigma nu Model --> egocentric noise only model
config.ModelName        =   "sigma_nu";
config.ParamName        =   ["sigma", "nu"];
config.NumParams        =   length(config.ParamName);
% Run the model
VAM

sigma_nu_IC             =   HealthyControls.Results.IC;

%% 2, beta, sigma, nu Model --> only encoding errors in distance
config.ModelName        =   "beta_sigma_nu";
config.ParamName        =   ["beta", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);
% Run the model
VAM

beta_sigma_nu_IC        =   HealthyControls.Results.IC;

%% 3, g2, g3, sigma, nu Model --> encoding errors and production errors in angle 
config.ModelName        =   "g2_g3_sigma_nu";
config.ParamName        =   ["g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);
% Run the model
VAM

g2_g3_sigma_nu_IC       =   HealthyControls.Results.IC;

%% 4, beta, g2, sigma, nu Model --> encoding error in distance and angle
config.ModelName        =   "beta_g2_sigma_nu";
config.ParamName        =   ["beta", "g2", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);
% Run the model
VAM

beta_g2_sigma_nu_IC     =   HealthyControls.Results.IC;

%% 5, beta, g3, sigma, nu Model --> encoding errors in distance and production error in angle
config.ModelName        =   "beta_g3_sigma_nu";
config.ParamName        =   ["beta", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);
% Run the model
VAM

beta_g3_sigma_nu_IC     =   HealthyControls.Results.IC;

%% 6, beta, g2, g3, sigma, nu Model --> encoding errors in distance and angles and production error in angle
config.ModelName        =   "beta_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);
% Run the model
VAM

beta_g2_g3_sigma_nu_IC  =   HealthyControls.Results.IC;

% %% 7, beta, g2, g3, k3, sigma, nu Model --> add another parameter k3 to the beta, g2, g3, sigma, nu Model
% config.ModelName        =   "beta_g2_g3_k3_sigma_nu";
% config.ParamName        =   ["beta", "g2", "g3", "k3", "sigma", "nu"];
% config.NumParams        =   length(config.ParamName);
% 
% Results_beta_g2_g3_k3_sigma_nu   =   getResultsAllConditions(YoungControls, config);
% beta_g2_g3_k3_sigma_nu_IC           =   Results_beta_g2_g3_k3_sigma_nu.IC;

%%
config.ResultFolder = pwd+"/Output/ModelFigures/ModelComparison_HC";;
%create storing folder for trajectory if not exist
if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

%%
ModelNames = {'M1', 'M2', 'M3', 'M4', 'M5', 'M6'};

%% Setting colors for using in plots
ColorPattern; 

%% reformating the IC structure
CondList = ["no change", "no distal cue", "no optical flow", "all"];
for idx_ = 1:4
    %cond = "no change"; 
    %cond = "no distal cue";
    %cond = "no optical flow";
    %cond = "all";
    cond = CondList(idx_);
    [sigma_nu_AIC,sigma_nu_BIC, sigma_nu_NLL] = reformatIC(sigma_nu_IC, cond);
    [beta_sigma_nu_AIC,beta_sigma_nu_BIC, beta_sigma_nu_NLL] = reformatIC(beta_sigma_nu_IC, cond);
    [g2_g3_sigma_nu_AIC,g2_g3_sigma_nu_BIC, g2_g3_sigma_nu_NLL] = reformatIC(g2_g3_sigma_nu_IC, cond);
    [beta_g2_sigma_nu_AIC, beta_g2_sigma_nu_BIC, beta_g2_sigma_nu_NLL] = reformatIC(beta_g2_sigma_nu_IC, cond);
    [beta_g3_sigma_nu_AIC,beta_g3_sigma_nu_BIC, beta_g3_sigma_nu_NLL] = reformatIC(beta_g3_sigma_nu_IC, cond);
    [beta_g2_g3_sigma_nu_AIC, beta_g2_g3_sigma_nu_BIC,beta_g2_g3_sigma_nu_NLL] = reformatIC(beta_g2_g3_sigma_nu_IC, cond);
    %[beta_g2_g3_k3_sigma_nu_AIC, beta_g2_g3_k3_sigma_nu_BIC,beta_g2_g3_k3_sigma_nu_NLL] = reformatIC(beta_g2_g3_k3_sigma_nu_IC, cond);
    
    % Box Plot of AIC
    ICType = "AIC";
    All_AIC = [sigma_nu_AIC, beta_sigma_nu_AIC, g2_g3_sigma_nu_AIC, beta_g2_sigma_nu_AIC, beta_g3_sigma_nu_AIC, beta_g2_g3_sigma_nu_AIC];
    plotBoxPlot(All_AIC, ModelNames, ICType, config, cond);
    
    % Box Plot of BIC
    ICType = "BIC";
    All_BIC = [sigma_nu_BIC, beta_sigma_nu_BIC, g2_g3_sigma_nu_BIC, beta_g2_sigma_nu_BIC, beta_g3_sigma_nu_BIC, beta_g2_g3_sigma_nu_BIC];
    plotBoxPlot(All_BIC, ModelNames, ICType, config, cond);
    
    % Box Plot of AIC
    ICType = "NegLogLikelihood";
    All_NLL = [sigma_nu_NLL, beta_sigma_nu_NLL, g2_g3_sigma_nu_NLL, beta_g2_sigma_nu_NLL, beta_g3_sigma_nu_NLL, beta_g2_g3_sigma_nu_NLL];
    plotBoxPlot(All_NLL, ModelNames, ICType, config, cond);

end
%% function for Box plot
function plotBoxPlot(data, ModelNames, ICType, config, cond)
    f = figure('visible','off','Position', [100 100 600 400]);
    
    %%%set paramsters
    num_boxplots = size(data,2);
    box_lineWidth = 0.5;
    whisker_value = 1.5;
    box_widths_value = 0.6;
    box_color_transparency = 0.5; %faceAlpha
    median_lineWidth = 2;
    median_color = 'k';
    scatter_jitter_value = 0.1;
    scatter_markerSize=10;
    scatter_marker_edgeColor = 'k';
    scatter_marker_edgeWidth = 0.5;
    scatter_color_transparency = 0.7; %faceAlpha

    %%% Font type and size setting %%%
    % Using Arial as default because all journals normally require the font to
    % be either Arial or Helvetica
    set(0,'DefaultAxesFontName','Arial')
    set(0,'DefaultTextFontName','Arial')
    set(0,'DefaultAxesFontSize',12)
    set(0,'DefaultTextFontSize',12)
   
    box_colors = config.color_scheme_npg(1:num_boxplots,:);
    
    %% main boxplot one box for each column in data
    bp = boxplot(data, 'whisker',whisker_value,'symbol','', ... %symbol ='' making outlier invisible
                'color','k', ...
                'Notch','on',...
                'widths',box_widths_value,...
                'labels', ModelNames);
    
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
        patch(get(h(i),'XData'),get(h(i),'YData'),box_colors(i,:),'FaceAlpha',box_color_transparency);
    end
    % Sending patch to back of the figure so that median can be drawn on top of it
    set(gca,'children',flipud(get(gca,'children'))) 

    %% Adjusting median
    h=findobj(gca,'tag','Median');
    for i = 1:length(h)
        h(i).LineWidth = median_lineWidth;
        h(i).Color = median_color;
    end

    %% add scatter plot and the mean of the data
    num_points = size(data,1);
    num_boxplots = size(data,2);
    for i=1:size(data,2)
        hold on
        x = i*ones(num_points,1)+scatter_jitter_value*(rand(num_points,1)-0.5); %jitter x
        scatter(x, data(:,i), scatter_markerSize, ...
                'filled','MarkerEdgeColor',scatter_marker_edgeColor, ...
                'MarkerFaceColor',box_colors(num_boxplots-i+1,:), ...
                'MarkerFaceAlpha',scatter_color_transparency,...
                'LineWidth',scatter_marker_edgeWidth); 
        hold on
        scatter(i, mean(data(:,i)), 4*scatter_markerSize, 'd',...
                'filled','MarkerEdgeColor','k', ...
                'MarkerFaceColor','w', ...
                'LineWidth',scatter_marker_edgeWidth);
    end

    %% Further post-processing the figure
%     if ICType=='AIC' | ICType=='BIC'
%         ylimup = 75;
%     else
%         ylimup = 35;
%     end    
    set(gca, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XColor'      , [.1 .1 .1], ...
        'YColor'      , [.1 .1 .1], ...
        'XTick'       , 1:1:100,... 
        'XLim'        , [0.5,length(ModelNames)+0.5],...
        'LineWidth'   , .5        );
    %'YLim'        , [0,ylimup],...
    title(cond)
    
    %% save figure
    exportgraphics(f,config.ResultFolder+"/"+cond+"_"+ICType+".png",'Resolution',300);
    exportgraphics(f,config.ResultFolder+"/"+cond+"_"+ICType+".pdf",'Resolution',300,'ContentType','vector');
end