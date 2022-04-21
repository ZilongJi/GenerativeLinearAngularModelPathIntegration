function Results = PerformGroupFit(GroupData, config)
%%Fit parameters for a single group
% GroupData is the extracted data containing the tracked positions of
% the participants as well as the their responses
% nsamples is the number of random samples to be taken when modelling the
% final position


TRIAL_FILTER = config.TrialFilter;          %load configurations necessary for the script

sampleSize = size(GroupData.FlagPos,2);     % Calculating sample size

% Initialize empty cell fro storing data
X               =       cell(1, sampleSize);% Actual positions
DX              =       cell(1, sampleSize);% Distances between subsequent points - segments li
THETADX         =       cell(1, sampleSize);% These are the angles between two subsequent segments. The angle indicate the rotation from the first segment towards the second, so it s the outer angle of the triangle.
segments        =       cell(1, length(sampleSize));
ProjSpeedL1     =       cell(1, sampleSize);% projected speed within detected start-to-end time window at leg 1
ProjSpeedL2     =       cell(1, sampleSize);% projected speed within detected start-to-end time window at leg 2
L1Dur           =       cell(1, sampleSize);% walking duration at leg 1
StandingDur     =       cell(1, sampleSize);% standing duration at cone2
L2Dur           =       cell(1, sampleSize);% walking duration at leg2
flagpos         =       cell(1, sampleSize);% flagpos
flagOoB         =       cell(1, sampleSize);% OoB flag
OoBLen          =       cell(1, sampleSize);% OoB length
GroupParameters =       cell(1, sampleSize);% Output value
IC              =       cell(1, sampleSize);  

anglebetween = @(va,vb) atan2d(va(:,1).*vb(:,2) - va(:,2).*vb(:,1), va(:,1).*vb(:,1) + va(:,2).*vb(:,2));

for j = 1:sampleSize
    if(TRIAL_FILTER == 0)
        %processing the data from all conditions
        flagpos{j}  = GroupData.FlagPos{j};
        OoBLen{j} = GroupData.Errors{j}.OoBLength;
        for idx = 1:size(GroupData.TrigPos{j},1)
            %If not out of bound or out of bound data is not present then take the trigpos
            if(GroupData.CondTable{j}.OutOfBound(idx) == 0 | isnan(GroupData.OutOfBoundPos{1,j}{idx}(1,1)))
                finalpos{j,1}(idx,1) = GroupData.TrigPos{j}(idx);
                flagOoB{j}(idx) = 0; %OoB flag is 0, i.e., not OoB trial
            elseif(GroupData.CondTable{j}.OutOfBound(idx) == 1)
                finalpos{j,1}(idx,1) = GroupData.ReconstructedOOB{j}.ReconstructedOoB(idx);
                flagOoB{j}(idx) = 1; %OoB flag is 0, i.e., it is OoB trial
            end
        end
    else
        %processing data according to "TRIAL_FILTER", with 1 no change, 2 no distal cues, 3 no optic flow
        flagpos{j}  = GroupData.FlagPos{j}(GroupData.CondTable{1,j}.Condition == TRIAL_FILTER);
        OoBLen{j}   = GroupData.Errors{j}.OoBLength(GroupData.CondTable{1,j}.Condition == TRIAL_FILTER);
        tempCnt     = 1;
        for idx = 1:size(GroupData.TrigPos{j},1)
            if(GroupData.CondTable{1,j}.Condition(idx) == TRIAL_FILTER)
                %If not out of bound or out of bound data is not present then take the trigpos
                if(GroupData.CondTable{j}.OutOfBound(idx) == 0 | isnan(GroupData.OutOfBoundPos{1,j}{idx}(1,1)))
                    finalpos{j,1}(tempCnt,1) = GroupData.TrigPos{j}(idx);
                    flagOoB{j}(tempCnt) = 0; %OoB flag is 0, i.e., not OoB trial
                elseif(GroupData.CondTable{j}.OutOfBound(idx) == 1)
                    finalpos{j,1}(tempCnt,1) = GroupData.ReconstructedOOB{j}.ReconstructedOoB(idx);
                    flagOoB{j}(tempCnt) = 1; %OoB flag is 1, i.e., it is OoB trial
                end
                tempCnt = tempCnt + 1;
            end
        end
    end

    if length(flagpos{j}) < config.NumFreeParams
        disp("%%%%%%%%%%%%%%% Skipping participant " + num2str(j) + ...
            ", because only "+ length(flagpos{j}) + ...
            " datapoints available for parameter estimation%%%%%%%%%%%%%%%\n");
        % NUmber of parameters in return matrix is fixed
        GroupParameters{j}  =   NaN(config.NumParams,1);
        IC{j}               =   nan;
        flagOoB{j}          =   [];
        continue;
    end

    % Create structures to save the data for each trial
    DX{j}           =   cell(1,length(flagpos{j}));
    THETADX{j}      =   cell(1,length(flagpos{j}));
    X{j}            =   cell(1,length(flagpos{j}));
    segments{j}     =   cell(1,length(flagpos{j}));
    ProjSpeedL1{j}  =   cell(2,length(flagpos{j}));
    ProjSpeedL2{j}  =   cell(2,length(flagpos{j}));
    L1Dur{j}        =   cell(1, length(flagpos{j}));
    StandingDur{j}  =   cell(1, length(flagpos{j}));
    L2Dur{j}        =   cell(1, length(flagpos{j}));
    
    for tr = 1:length(flagpos{j})
        
        X{j}{tr}         =      [flagpos{j}{tr}(:,[1,3]);finalpos{j}{tr}([1,3])];
        segments{j}{tr}  =      X{j}{tr}(2:end,:) - X{j}{tr}(1:end-1,:);
        DX{j}{tr}        =      sqrt(sum(segments{j}{tr}.^2,2));        
     
        outer_rad        =      deg2rad([0; anglebetween(segments{j}{tr}(1:end-1,:), segments{j}{tr}(2:end,:))]);   
        THETADX{j}{tr}   =      mod(outer_rad, 2*pi);%wrap the angle into (0,2pi)

        %extract the projected speed information along with the time information on outbound path
        L1_Vel_proj             =       GroupData.TrackedL1{j}{tr}.Vel_proj;
        L1_Time                 =       GroupData.TrackedL1{j}{tr}.Time;
        L1_Filtered_Vel_proj    =       GroupData.TrackedL1{j}{tr}.Filtered_Vel_proj;
        L1_Vel_proj_selected    =       L1_Vel_proj(L1_Filtered_Vel_proj);
        L1_Time_selected        =       L1_Time(L1_Filtered_Vel_proj);
        ProjSpeedL1{j}{1,tr}    =       L1_Time_selected;  %selected time in a start-to-end range
        ProjSpeedL1{j}{2,tr}    =       L1_Vel_proj_selected;  %selected speed in a start-to-end range

        L2_Vel_proj             =       GroupData.TrackedL2{j}{tr}.Vel_proj;
        L2_Time                 =       GroupData.TrackedL2{j}{tr}.Time;
        L2_Filtered_Vel_proj    =       GroupData.TrackedL2{j}{tr}.Filtered_Vel_proj;
        L2_Vel_proj_selected    =       L2_Vel_proj(L2_Filtered_Vel_proj);
        L2_Time_selected        =       L2_Time(L2_Filtered_Vel_proj);
        ProjSpeedL2{j}{1, tr}   =       L2_Time_selected;  %selected time in a start-to-end range
        ProjSpeedL2{j}{2, tr}   =       L2_Vel_proj_selected; %selected speed in a start-to-end range

        %extract the standing duration at cone2
        Idx_Cond                =       GroupData.CondTable{j}.Condition == TRIAL_FILTER;
        standing_duration       =       GroupData.Reconstructed{j}.T_Standing(Idx_Cond);
        StandingDur{j}{tr}      =       standing_duration(tr);
    end


    %put all of things we need into a struct for sending to FitData
    Input.DX            =   DX{j};
    Input.THETADX       =   THETADX{j};
    Input.X             =   X{j};
    Input.ProjSpeedL1   =   ProjSpeedL1{j};
    Input.ProjSpeedL2   =   ProjSpeedL2{j};
    Input.StandingDur   =   StandingDur{j};

    % Do the data fitting
    disp(['%%%%%%%%%%%%%%% STARTING FIT PER PARTICIPANT ' num2str(j) ' %%%%%%%%%%%%%%%']);
    [GroupParameters{j}, IC{j}] = FitData(DX{j}, THETADX{j}, X{j}, ProjSpeedL1{j}, ProjSpeedL2{j}, config);
end

%Transforming the fitted parameters to array
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

end