function [negloglikelihood] = Estimate_beta_k_g2_sigma_nu(beta, k, g2, sigma, nu, Input, config)
% Zilong Ji, UCL, 2022, zilong.ji@ucl.ac.uk
% Find the likelihood with parameters specified in the model function name
% Args (in common between all estimates functions):
% beta: decay factor for the mental distance (calculation error)
% k: velocity gain factor for the leaky integrator (encoding error)
% g2: rotation gain for the second turn (encoding error)
% g3: regression to the mean effect in return angle (production error)
% m3: regression to mean effect in return distance (production error)
% sigma: standard deviation for the Gaussian distribution of the return
% distance (calculation + unexplained error)
% nu: standard deviation for the Gaussian distribution of the return angle
% ((calculation + unexplained error))
% Input: contains all the data information for estimating, see PerformGroupFit for how it was generated
% For schematics about the steps of the generative linear and angular
% model please refer to Fig. 2 and online methods
% conif: please refer to GLAMPI_PrepareBaseConfig
% ===================================================================================

%% information necessary for running parameter estimation
DX              =   Input.DX;
THETAX          =   Input.THETADX;
L1Dur           =   Input.L1Dur;
L2Dur           =   Input.L2Dur;
StandingDur     =   Input.StandingDur;
flagOoB         =   Input.flagOoB;

sampleSize          =   size(DX,2);
negloglikelihood    =   0;

for tr = 1:sampleSize
    %% extract the physical data info
    l1          =       DX{tr}(1);
    l2          =       DX{tr}(2);
    l3          =       DX{tr}(3);
    theta2      =       THETAX{tr}(2); 
    theta3      =       THETAX{tr}(3); 
    durationL1  =       L1Dur{tr}; 
    durationL2  =       L2Dur{tr};
    durationStand =     StandingDur{tr};
    
    %% whether to use weber's law to scaling the noise strength
    if config.useweber == true
        scale = sqrt(l1^2+l2^2); %noise scale
    else
        scale = 1; %noise scale
    end
    sigma_scaled = scale*sigma;
    nu_scaled = scale*nu;

    %mental point 1 (asuming a constant speed)
    %considering standing duration or not
    if config.includeStand==true
        men_length1 = l1*k*(1-exp(-beta*durationL1))/(beta*durationL1)*exp(-beta*(durationL2+durationStand));
    else
        men_length1 = l1*k*(1-exp(-beta*durationL1))/(beta*durationL1)*exp(-beta*durationL2);
    end
    men_p1 = [men_length1,0];
    
    theta2_prime = g2*theta2;

    %mental point 2, (asuming a constant speed)
    men_length2 = l2*k*(1-exp(-beta*durationL2))/(beta*durationL2);
    men_p2      = [men_length1+men_length2*cos(theta2_prime),men_length2*sin(theta2_prime)];

    %calculate length of mental vector 3
    h           = norm(men_p2);
    %calculate turn angle of mental vector 3
    vec1        = men_p2-men_p1; 
    vec2        = [0,0]-men_p2;
    alpha       = atan2d(vec1(1)*vec2(2)-vec1(2)*vec2(1),vec1(1)*vec2(1)+vec1(2)*vec2(2));
    alpha       = deg2rad(alpha);   %transfer from degree to radians
    
    %mental return angle
    theta3_prime = alpha;
    
    %angular noise difference
    angluar_diff = theta3-theta3_prime;
    %the negative loglikelihood of angle
    %neg_ll_angle = log(2*pi) + log(besseli(0,nu_scaled)) - nu_scaled*cos(angluar_diff);
    neg_ll_angle = 1/2*log(2*pi) + log(nu_scaled) + (angluar_diff^2)/(2*nu_scaled^2);

    %distance noise difference
    l3_prime    = h;
    dist_diff   = l3-l3_prime;

    %the negative loglikelihood of distance on non-OoB trials
    if flagOoB(tr)==0
        %this is a non-OoB trial
        neg_ll_dist = 1/2*log(2*pi) + log(sigma_scaled) + (dist_diff^2)/(2*sigma_scaled^2);
    else
        %this is an OoB trial
        neg_ll_dist = 0;
    end

    %total negative loglikelihood
    neg_ll = neg_ll_angle + neg_ll_dist;

    negloglikelihood = negloglikelihood + neg_ll;
end
end