function [kValue_lh, concentration_lh, lhSample] = latinHypercuber_240906(kValue_vector, concentration_vector, trials, ...
    kValue_range, concentration_range, kValue_noChange, concentration_noChange)
    % Function to create a random parameter set for the MPSA using the latin
    % hypercube sampling. The k-values and concentrations are concatenated into
    % a single vector in the function for parameter sampling. Following
    % parameter sampling, the k-value and concentration vectors are outputted
    % as individual vectors for the MPSA run. 
    %   kValue_lh = latin hypercube sampled k-value parameters that are
    %   inputted in the MPSA
    %   concentration_lh = latin hypercube sampled concentration parameters
    %   that are inputted in the MPSA    
    %   kValue_vector = calibrated k-values
    %   concentration_vector = calibrated initial concentrations
    %   trials = number of iterations for the MPSA to run
    %   kValue_range = non-log scale fold change to create range to sample parameters
    %   concentration_range = non-log scale fold change to create range to sample parameters
    %   kValue_noChange = k-values that are not changed in the parameter
    %   sampling
    %   concentration_noChange = concentrations that are not changed in the
    %   parameter sampling

    % 240908 update: allowing for specific ranges for parameters, i.e., not
    % all k-values or concentrations assume the same parameter range in
    % the MPSA

    
% combining k-values and concentrations into a single array (calibrated
% values)
combinedParameters = [kValue_vector', concentration_vector'];

% % Range for LH sampling
% kValue_range = [10^-log10(kValue_range), 10^log10(kValue_range)];
%     kValue_rangeLow = kValue_range(1);
%     kValue_rangeHigh = kValue_range(2);
% conc_range = [10^-log10(concentration_range), 10^log10(concentration_range)];
%     conc_rangeLow = conc_range(1);
%     conc_rangeHigh = conc_range(2);

% Specific range for each parameter
kValue_range_max = kValue_range(1,:); % max FC
kValue_range_min = kValue_range(2,:); % min FC
concentration_range_max = concentration_range(1,:); % max FC
concentration_range_min = concentration_range(2,:); % min FC

% kValue_range; % includes the max and min FC
% concentration_range; % includes the max and min FC

% latin hypercube sampling (0,1)
lhSample = lhsdesign(trials,length(combinedParameters)); % number of replications (rows), number of parameters (columns)

% factoring in LH sampling to the range of k-value and concentration values
% -> shifting the mean to 0 instead of 0.5
combinedParameters_logFactor = nan(size(lhSample));

for k = 1:length(combinedParameters)

    % locating columns in the combined array that represent k-values;
    % calculate log-scale fold change for each parameter for each iteration
    if k <= length(kValue_vector)
%         combinedParameters_logFactor(:,k) = lhSample(:,k).*log10(kValue_rangeHigh./kValue_rangeLow)+log10(kValue_rangeLow);
        combinedParameters_logFactor(:,k) = lhSample(:,k).*log10(kValue_range_max(k)./kValue_range_min(k))+log10(kValue_range_min(k));

    % locating columns in the combined array that represent concentrations;
    % calculate log-scale fold change for each parameter for each iteration
    elseif k > length(kValue_range)
%         combinedParameters_logFactor(:,k) = lhSample(:,k).*log10(conc_rangeHigh./conc_rangeLow)+log10(conc_rangeLow);
        % concentratoin_range_... vector is of length 40. k >
        % length(kValue_range) (i.e., number of kValues) starts looking at
        % concentrations. I need to add length(kValue_range) to each
        % concentration "k" to select the correct value
        k_concCorrect = k - length(kValue_range);
        combinedParameters_logFactor(:,k) = lhSample(:,k).*log10(concentration_range_max(k_concCorrect)./concentration_range_min(k_concCorrect))+log10(concentration_range_min(k_concCorrect)); 

    end
end

% Convert log scale factor to linear scale
combinedParameters_linearScale = 10.^combinedParameters_logFactor;

% factor in parameter values to create parameter sets for MPSA
combinedParameters_mpsa = combinedParameters_linearScale.*combinedParameters;

% replacing parameter value to their calibrated value if selected to not 
% change in the MPSA
species_noChange = [kValue_noChange, (concentration_noChange+length(kValue_vector))]; %combined non-changed kValue and concentration parameters into a single list

for n = 1:length(species_noChange)
    parameterNum = species_noChange(n);
    combinedParameters_mpsa(:,parameterNum) = combinedParameters(parameterNum);
end

% outputting parameter values in k-value and concentration specific arrays
kValue_lh = combinedParameters_mpsa(:,1:length(kValue_vector));
kValue_lh = kValue_lh';

concentration_lh = combinedParameters_mpsa(:,1+length(kValue_vector):length(kValue_vector)+length(concentration_vector));
concentration_lh = concentration_lh';

end




