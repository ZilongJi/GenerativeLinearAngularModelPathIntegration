%% Cleaning variables and set intial seed for code reproducibility
clearvars; close all; clc;
rng('default'); 

%% Loading data
disp('%%%%%%%%%%%%%%% DATA LOADING ... %%%%%%%%%%%%%%%');
load('AllDataErrorsPrevent.mat');

%% setting the configuration
config.Speed.alpha                                      = 0.9;    % Paramanter for running speed calculation
config.Speed.timeOffsetAfterFlagReach                   = 0.5;    % Time to track after cone reached in seconds 
config.Speed.smoothWindow                               = 20;     % tracking rate should be 10Hz so 4 secs window is 40 datapoints
config.Speed.velocityCutoff                             = 0.2;    % velocity cutoff to select only the walking part of the reconstructed velocity
config.Speed.timeOffsetForDetectedTemporalWindow        = 0.2;    % time in seconds that will push earlier/ the detected rising edge
config.TrackedInboundAngularDeltaT                      = 1;      % delta time step to integrate the angular information from the tracking data
config.includeStand                                     = false;  % set to true if the model should also include for duration of time spent on each segment the time calculate from tracking data while standing at each cone
config.useOoBtrials                                     = true;   % Are we using out of bound trials in the model?
config.useTrialFilter                                   = true;   % when true the model will be fitted for each of the task conditions separately. If false it will discard the
config.Speed.tresholdForBadParticipantL1Recontruction   = 0.0;    % threshold for escluding participants with the weird shaped trials (on l1). If zero all data will be used.
disp('%%%%%%%%%%%%%%% Data loading complete... %%%%%%%%%%%%%%%');

%% Preprocessing FHPos
disp("%%%%%%%%%%%%%%% Preprocessing FHPos %%%%%%%%%%%%%%%");
FamilyHistPos   = TransformPaths(FamilyHistPos);%transform data
FamilyHistPos   = CalculateTrackingPath(FamilyHistPos, config);
ManuallyScoringFamilyHistPos;

%% Preprocessing FHNeg
disp("%%%%%%%%%%%%%%% Preprocessing FHNeg %%%%%%%%%%%%%%%");
FamilyHistNeg   = TransformPaths(FamilyHistNeg);%transform data
FamilyHistNeg   = CalculateTrackingPath(FamilyHistNeg, config);
ManuallyScoringFamilyHistNeg;

% Model fitting
config.ModelName        =   "beta_k_g2_g3_sigma_nu";
config.ParamName        =   ["beta", "k", "g2", "g3", "sigma", "nu"];
config.NumParams        =   length(config.ParamName);

FamilyHistPos.Results = getResultsAllConditions(FamilyHistPos, config);
FamilyHistNeg.Results = getResultsAllConditions(FamilyHistNeg, config);

%% Preparing the output
config.ResultFolder     =   pwd + "/Output/Coco/"+config.ModelName+"/FourGroup";

if ~exist(config.ResultFolder, 'dir')
   mkdir(config.ResultFolder);
end

% Generating color scheme for our paper
ColorPattern; 

% Collecting information from output
AllFamilyHistPosParams     =   FamilyHistPos.Results.estimatedParams;
AllFamilyHistNegParams     =   FamilyHistNeg.Results.estimatedParams;


%% 4 groups FH+ Apoe+ v.s. FH+ Apoe- v.s. FH- Apoe+ v.s. FH- Apoe-

%find the apoe tag for FH+ and FH-
FHPos_Apoe = zeros(length(FamilyHistPos.Info),1);
for i=1:length(FamilyHistPos.Info)
    name = FamilyHistPos.Info{i};
    idx = find(strcmp(preventTab.subjectID, name));
    %find the apoe tag
    apoe_label = preventTab.apoe4(idx);
    FHPos_Apoe(i) = apoe_label;
end

FHNeg_Apoe = zeros(length(FamilyHistNeg.Info),1);
for i=1:length(FamilyHistNeg.Info)
    name = FamilyHistNeg.Info{i};
    idx = find(strcmp(preventTab.subjectID, name));
    %find the apoe tag
    apoe_label = preventTab.apoe4(idx);
    FHNeg_Apoe(i) = apoe_label;
end

%extract the subgroups 1,FH+ APoe+ 2,FH+ APoe-....
Condition = [1,2,3];  %switch the no distal cue condition with the non optical flow condition
for k=1:3
    cond = Condition(k);
    Params = AllFamilyHistPosParams{1,cond};
    Gender = FamilyHistPos.Gender;

    nonNanIdx = ~isnan(FHPos_Apoe);
    nonNanParams = Params(nonNanIdx,:);
    nonNanApoe = FHPos_Apoe(nonNanIdx);
    nonNanGender = Gender(nonNanIdx);

    FHPos_ApoePos_Param{1,k} = nonNanParams(logical(nonNanApoe),:);
    FHPos_ApoeNeg_Param{1,k} = nonNanParams(~logical(nonNanApoe),:);
    FHPos_ApoePos_Gender = nonNanGender(logical(nonNanApoe));
    FHPos_ApoeNeg_Gender = nonNanGender(~logical(nonNanApoe));
end

%extract the subgroups 1,FH+ APoe+ 2,FH+ APoe-....
for k=1:3
    cond = Condition(k);
    Params = AllFamilyHistNegParams{1,cond};
    Gender = FamilyHistNeg.Gender;

    nonNanIdx = ~isnan(FHNeg_Apoe);
    nonNanParams = Params(nonNanIdx,:);
    nonNanApoe = FHNeg_Apoe(nonNanIdx);
    nonNanGender = Gender(nonNanIdx);

    FHNeg_ApoePos_Param{1,k} = nonNanParams(logical(nonNanApoe),:);
    FHNeg_ApoeNeg_Param{1,k} = nonNanParams(~logical(nonNanApoe),:);
    FHNeg_ApoePos_Gender = nonNanGender(logical(nonNanApoe));
    FHNeg_ApoeNeg_Gender = nonNanGender(~logical(nonNanApoe));
end

config.FHPos_ApoePos_Gender = FHPos_ApoePos_Gender;
config.FHPos_ApoeNeg_Gender = FHPos_ApoeNeg_Gender;
config.FHNeg_ApoePos_Gender = FHNeg_ApoePos_Gender;
config.FHNeg_ApoeNeg_Gender = FHNeg_ApoeNeg_Gender;

% three way anova on gender, codition, and 4 FH+APOE groups
[anova_tab,...
 multicomp_tab1,...
 multicomp_tab2, ...
 multicomp_tab3, ...
 multicomp_tab12, ...
 multicomp_tab13, ...
 multicomp_tab23, ...
 multicomp_tab123] = ThreewayAnova_4group(FHPos_ApoePos_Param, ...
                                                           FHPos_ApoeNeg_Param, ...
                                                           FHNeg_ApoePos_Param, ...
                                                           FHNeg_ApoeNeg_Param, ...
                                                           config);

%% two way anova
[anova_tab,multicomp_tab1,~, ~] = TwowayAnova_4group(FHPos_ApoePos_Param, FHPos_ApoeNeg_Param, FHNeg_ApoePos_Param, FHNeg_ApoeNeg_Param,config);