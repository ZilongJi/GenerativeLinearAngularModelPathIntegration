function Results = PerformGroupFit(GroupData, config)
%% PerformGroupFit
% Andrea Castegnaro, UCL, uceeaca@ucl.ac.uk
% Zilong Ji, UCL, zilong.ji@ucl.ac.uk
% Fit parameters for a single group. Calculates behavioural results. 
% GroupData is the data obtained at the end of the pre-processing stage.
% ===================================================================================

TRIAL_FILTER = config.TrialFilter;          % check wheter run the script by keeping environmental conditions separate or not
subjectNum = size(GroupData.FlagPos,2);     

%% Initialize empty cell for storing output data
X               =       cell(1, subjectNum); % Actual cone positions
DX              =       cell(1, subjectNum); % Distances between subsequent locations - e.g. segments li for outbound and inbound path
THETADX         =       cell(1, subjectNum); % Angles between two subsequent locations. Always the outer angle of the triangle (see Fig1 for a schematic). For the last angle we calculated the egocentric rotation of the participant using the tracking data
segments        =       cell(1, subjectNum); % Vector differences between subsequent locations (between cones and between triggered position and last cone).
correctReDist   =       cell(1, subjectNum); % Length of vector cone 3 - 1 - used only locally
correctReAngle  =       cell(1, subjectNum); % Outer Angle between vector cone 2-3 and vector cone 3-1 (short angle)
DistErr         =       cell(1, subjectNum); % Difference between correct return distance and walked participant distance (length vector triggered location - cone 3)
AngleErr        =       cell(1, subjectNum); % Difference between correct return angle and participant rotation at cone 3
PropDistErr     =       cell(1, subjectNum); % Ratio between length vector cone 3 - triggered position and length vector cone 3 - cone 1
PropAngErr      =       cell(1, subjectNum); % Ratio between participant s rotation at cone 3 and real returning angle
LocationErr     =       cell(1, subjectNum); % Length of vector cone 3 - triggered location
ProjSpeedL1     =       cell(1, subjectNum); % Projected speed within detected start-to-end time window at leg 1
ProjSpeedL2     =       cell(1, subjectNum); % Projected speed within detected start-to-end time window at leg 2
L1Dur           =       cell(1, subjectNum); % Walking duration at leg 1
L2Dur           =       cell(1, subjectNum); % Walking duration at leg2
StandingDur     =       cell(1, subjectNum); % Standing duration at cone 2
flagpos         =       cell(1, subjectNum); % Exctracting flag positions from input data
flagOoB         =       cell(1, subjectNum); % Flagging OoB trials
GroupParameters =       cell(1, subjectNum); % The values of the fitted parameters for a given model
IC              =       cell(1, subjectNum); % Information Criterion (aic, bic and (neg) likelihood) after fitting the model
TRIALs          =       zeros(subjectNum, 3); % Number of trials

%% help function to calculate the angle between two vector
anglebetween = @(va,vb) atan2d(va(:,1).*vb(:,2) - va(:,2).*vb(:,1), va(:,1).*vb(:,1) + va(:,2).*vb(:,2));

% Looping through each of the participant
for j = 1:subjectNum

    % Filter out participants who did short inbound walking, or that
    % retraced the triangle to get back to cone 1, or that did not have out
    % of bound location registered from the tracking data. See
    % preprocessing for details and online methods section.

    if ismember(j, GroupData.BadPptIdxs)
        % set results 
        GroupParameters{j}  =   NaN(config.NumParams,1);
        IC{j}.aic           =   nan;
        IC{j}.bic           =   nan;
        IC{j}.negll         =   nan;
        IC{j}.likelihood    =   nan;
        flagOoB{j}          =   [];
        disp(['%%%%%%%%%%%%%%% Skipping PARTICIPANT ' num2str(j) ' ---- because excluded after preprocessing%%%%%%%%%%%%%%%']);
        continue
    end
    
    if(TRIAL_FILTER == 0)
        % We don t fit based on the environmental condition
        % Merging all three environmental conditions
        flagpos{j}  = GroupData.FlagPos{j};
        BadExecutionTrials  = GroupData.Reconstructed{j}.BadExecution;
        realReturnAngles    = GroupData.Reconstructed{j}.RealReturnAngle;
        TrialNum    = size(GroupData.TrigPos{j},1);
        for idx = 1:TrialNum
            % Getting final position based on the fact this is an out of
            % bound trial or not (please see CalculateOoB for details)
            if(GroupData.CondTable{j}.OutOfBound(idx) == 0 | isnan(GroupData.OutOfBoundPos{1,j}{idx}(1,1)))
                finalpos{j}(idx) = GroupData.TrigPos{j}(idx);
                flagOoB{j}(idx) = 0; % not OoB trial
            elseif(GroupData.CondTable{j}.OutOfBound(idx) == 1)
                finalpos{j}(idx) = GroupData.ReconstructedOOB{j}.ReconstructedOoB(idx);
                flagOoB{j}(idx) = 1; % OoB trial
            end
        end   
        Idx_GoodTrials      = BadExecutionTrials == 0;
        flagpos{j}          = flagpos{j}(Idx_GoodTrials);
        realReturnAngles    = realReturnAngles(Idx_GoodTrials);
        finalpos{j}         = finalpos{j}(Idx_GoodTrials);
        flagOoB{j}          = flagOoB{j}(Idx_GoodTrials); 
    
        leg1_duration       = GroupData.Reconstructed{j}.T_L1(Idx_GoodTrials);          %filter the duration of subject j at leg 1
        leg2_duration       = GroupData.Reconstructed{j}.T_L2(Idx_GoodTrials);          %filter the duration of subject j at leg 2
        standing_duration   = GroupData.Reconstructed{j}.T_Standing(Idx_GoodTrials);    %filter the duration of subject j standing at cone2 
        TrackedL1           = GroupData.TrackedL1{j}(Idx_GoodTrials);
        TrackedL2           = GroupData.TrackedL2{j}(Idx_GoodTrials);
    else
        %%filter trails based on the TRIAL_FILTER, with 1 no change, 2 no distal cues, 3 no optic flow
        Idx_Cond            = GroupData.CondTable{j}.Condition == TRIAL_FILTER; 
        flagpos{j}          = GroupData.FlagPos{j}(Idx_Cond);
        BadExecutionTrials  = GroupData.Reconstructed{j}.BadExecution(Idx_Cond);
        realReturnAngles    = GroupData.Reconstructed{j}.RealReturnAngle(Idx_Cond);
    
        tempCnt             = 1;
        TrialNum            = size(GroupData.TrigPos{j},1);
    
        %get the final position and OoB flag 
        for idx = 1:TrialNum %for each trial, if belongs to TRIAL_FILTER, go into
            if(GroupData.CondTable{j}.Condition(idx) == TRIAL_FILTER)
                %If not out of bound or out of bound data is not present then take the trigpos
                if(GroupData.CondTable{j}.OutOfBound(idx) == 0 | isnan(GroupData.OutOfBoundPos{1,j}{idx}(1,1)))
                    finalpos{j}(tempCnt) = GroupData.TrigPos{j}(idx);
                    flagOoB{j}(tempCnt) = 0; %OoB flag is 0, i.e., not OoB trial
                elseif(GroupData.CondTable{j}.OutOfBound(idx) == 1)
                    finalpos{j}(tempCnt) = GroupData.ReconstructedOOB{j}.ReconstructedOoB(idx);
                    flagOoB{j}(tempCnt) = 1; %OoB flag is 1, i.e., it is OoB trial
                end
                tempCnt = tempCnt + 1;
            end
        end
        
        Idx_GoodTrials      = BadExecutionTrials == 0;
        flagpos{j}          = flagpos{j}(Idx_GoodTrials);
        realReturnAngles    = realReturnAngles(Idx_GoodTrials);
        finalpos{j}         = finalpos{j}(Idx_GoodTrials);
        flagOoB{j}          = flagOoB{j}(Idx_GoodTrials); 
    
        leg1_duration       = GroupData.Reconstructed{j}.T_L1(Idx_Cond);          %filter the duration of subject j at leg 1
        leg1_duration       = leg1_duration(Idx_GoodTrials);  
    
        leg2_duration       = GroupData.Reconstructed{j}.T_L2(Idx_Cond);          %filter the duration of subject j at leg 2
        leg2_duration       = leg2_duration(Idx_GoodTrials);
    
        standing_duration   = GroupData.Reconstructed{j}.T_Standing(Idx_Cond);    %filter the duration of subject j standing at cone2 
        standing_duration   = standing_duration(Idx_GoodTrials);
    
        TrackedL1           = GroupData.TrackedL1{j}(Idx_Cond);
        TrackedL1           = TrackedL1(Idx_GoodTrials);
    
        TrackedL2           = GroupData.TrackedL2{j}(Idx_Cond);
        TrackedL2           = TrackedL2(Idx_GoodTrials);
    end
    
    %get X, segments, DX etc....
    for tr = 1:length(flagpos{j})
        
        X{j}{tr}         =      [flagpos{j}{tr}(:,[1,3]);finalpos{j}{tr}([1,3])];
        segments{j}{tr}  =      X{j}{tr}(2:end,:) - X{j}{tr}(1:end-1,:);
        % DX contains the lenght of each segment of the outbound path plus the
        % inbound path
        DX{j}{tr}        =      sqrt(sum(segments{j}{tr}.^2,2));        
     
        outer_rad        =      [0; anglebetween(segments{j}{tr}(1:end-1,:), segments{j}{tr}(2:end,:))];
        outer_rad(3,1)   =      realReturnAngles(tr);     
        THETADX{j}{tr}   =      deg2rad(outer_rad);
        
        if config.useOoBtrials == true
            
            % Calculating triangle property
            p1 = X{j}{tr}(1,:);
            p2 = X{j}{tr}(2,:);
            p3 = X{j}{tr}(3,:);
            tp = X{j}{tr}(4,:);

            vec1 = p3-p2; %second leg
            vec2 = p1-p3; %third leg (real)
            rvec = tp-p3; %participant return vector for current trial
            %extract correct return distance and angle
            %consider InB trials and the angular information in OoB trails
            if flagOoB{j}(tr) == 0
                x_2 = X{j}{tr}(3,:);
                correctReDist{j}{tr} = sqrt(sum(x_2.^2));
                DistErr{j}{tr}       = correctReDist{j}{tr}-DX{j}{tr}(3);
                realReturnLength     = norm(rvec);
                PropDistErr{j}{tr}   = realReturnLength/correctReDist{j}{tr};
                x_3 = X{j}{tr}(4,:);   %trigger position
                LocationErr{j}{tr}   = sqrt(sum(x_3.^2));
            else 
                DistErr{j}{tr}       = nan;
                PropDistErr{j}{tr}   = nan;
                LocationErr{j}{tr}   = nan;
            end    

            correctReAngle{j}{tr} = anglebetween(vec1, vec2); 
            
            AngleErr{j}{tr}   = wrapTo180(realReturnAngles(tr)-correctReAngle{j}{tr});
            PropAngErr{j}{tr} = wrapTo360(realReturnAngles(tr))/correctReAngle{j}{tr};
        else
            %use only InB trials
            %extract correct return distance and angle
            %consider only InB trials
            if flagOoB{j}(tr) == 0
                x_2 = X{j}{tr}(3,:);
                correctReDist{j}{tr} = sqrt(sum(x_2.^2));
    
                p1 = X{j}{tr}(1,:);
                p2 = X{j}{tr}(2,:);
                p3 = X{j}{tr}(3,:);
                tp = X{j}{tr}(4,:);
                vec1 = p3-p2; 
                vec2 = p1-p3;
                rvec = tp-p3;
                correctReAngle{j}{tr} = anglebetween(vec1, vec2);
                realReturnLength     = norm(rvec);
    
                %calculate distance error and angular error
                DistErr{j}{tr}     = correctReDist{j}{tr}-DX{j}{tr}(3);
                PropDistErr{j}{tr} = realReturnLength/correctReDist{j}{tr};
                %no wrap
                %AngleErr{j}{tr} = abs(correctReAngle{j}{tr}-realReturnAngles(tr));
                %wrap the real return angle to [0,360]
                %AngleErr{j}{tr} = correctReAngle{j}{tr}-wrapTo360(realReturnAngles(tr));
                %wrap the angular error to [-180,180]
                %AngleErr{j}{tr}   = wrapTo180(correctReAngle{j}{tr}-realReturnAngles(tr));
                AngleErr{j}{tr}   = wrapTo180(realReturnAngles(tr)-correctReAngle{j}{tr});
                PropAngErr{j}{tr} = wrapTo360(realReturnAngles(tr))/correctReAngle{j}{tr};

                x_3 = X{j}{tr}(4,:);   %trigger position
                LocationErr{j}{tr}   = sqrt(sum(x_3.^2));         
            else 
                DistErr{j}{tr}     = nan;
                PropDistErr{j}{tr} = nan;
                AngleErr{j}{tr}    = nan;
                PropAngErr{j}{tr}  = nan;
                LocationErr{j}{tr} = nan;
            end
        end

        %extract the projected speed information along with the time information on outbound path
        L1_Vel_proj             =       TrackedL1{tr}.Vel_proj;
        L1_Time                 =       TrackedL1{tr}.Time;
        L1_Filtered_Vel_proj    =       TrackedL1{tr}.Filtered_Vel_proj;
        L1_Vel_proj_selected    =       L1_Vel_proj(L1_Filtered_Vel_proj);
        L1_Time_selected        =       L1_Time(L1_Filtered_Vel_proj);
        ProjSpeedL1{j}{1,tr}    =       L1_Time_selected;  %selected time in a start-to-end range
        ProjSpeedL1{j}{2,tr}    =       L1_Vel_proj_selected;  %selected speed in a start-to-end range

        L2_Vel_proj             =       TrackedL2{tr}.Vel_proj;
        L2_Time                 =       TrackedL2{tr}.Time;
        L2_Filtered_Vel_proj    =       TrackedL2{tr}.Filtered_Vel_proj;
        L2_Vel_proj_selected    =       L2_Vel_proj(L2_Filtered_Vel_proj);
        L2_Time_selected        =       L2_Time(L2_Filtered_Vel_proj);
        ProjSpeedL2{j}{1, tr}   =       L2_Time_selected;  %selected time in a start-to-end range
        ProjSpeedL2{j}{2, tr}   =       L2_Vel_proj_selected; %selected speed in a start-to-end range
        
        L1Dur{j}{tr}            =       leg1_duration(tr);      %extract the walking duration at leg 1 
        L2Dur{j}{tr}            =       leg2_duration(tr);      %extract the walking duration at leg 2
        StandingDur{j}{tr}      =       standing_duration(tr);  %extract the standing duration at cone2
    end


    %% put all of the things we need into a struct for sending to FitData
    Input.DX                 =   DX{j};
    Input.THETADX            =   THETADX{j};
    Input.X                  =   X{j};
    Input.flagOoB            =   flagOoB{j};
    Input.ProjSpeedL1        =   ProjSpeedL1{j};
    Input.ProjSpeedL2        =   ProjSpeedL2{j};
    Input.L1Dur              =   L1Dur{j};
    Input.L2Dur              =   L2Dur{j};
    Input.StandingDur        =   StandingDur{j};
    
    %record the number of trials

    TRIALs(j, TRIAL_FILTER) = length(flagpos{j});
    if length(flagpos{j}) < config.NumParams-1
        %lack of trials, skip estimation
        disp("%%%%%%%%%%%%%%% Skipping participant " + num2str(j) + ...
            ", because only "+ length(flagpos{j}) + ...
            " datapoints available for parameter estimation%%%%%%%%%%%%%%%\n");
        % set results to nan for later processing
        GroupParameters{j}  =   NaN(config.NumParams,1);
        IC{j}.aic           =   nan;
        IC{j}.bic           =   nan;
        IC{j}.negll         =   nan;
        IC{j}.likelihood    =   nan;
        %flagOoB{j}          =   [];
        continue;
    end

    %% Do the data fitting
    disp(['%%%%%%%%%%%%%%% STARTING FIT PER PARTICIPANT ' num2str(j) ' %%%%%%%%%%%%%%%']);
    [GroupParameters{j}, IC{j}] = FitData(Input, config);
end

%%
%Transforming the fitted parameters from cell to array
[~, rows]       = size(GroupParameters);
[cols,~]        = size(GroupParameters{1});
estimatedParams = zeros(rows,cols);
for i=1:rows
   estimatedParams(i,:) = GroupParameters{i}';
end

%put all results into a matlab structure
Results.estimatedParams =   estimatedParams;
Results.X               =   X;
Results.DX              =   DX;
Results.THETADX         =   THETADX;
Results.IC              =   IC;
Results.flagOoB         =   flagOoB;
Results.DistErr         =   DistErr;
Results.AngleErr        =   AngleErr;
Results.PropDistErr     =   PropDistErr;
Results.PropAngErr      =   PropAngErr;
Results.LocationErr     =   LocationErr;
Results.L1Dur           =   L1Dur;
Results.L2Dur           =   L2Dur;
Results.TRIALs          =   TRIALs;

end