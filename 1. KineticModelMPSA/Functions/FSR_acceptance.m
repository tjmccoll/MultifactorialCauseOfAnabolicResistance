function [pass] = FSR_acceptance(dataFSR, FSR_range, t, x, eqDur, simDur_end, species)
    % Function assesses the acceptance of the FSR integral of the MPSA 
    % simulation
    %   dataFSR = the FSR integral calculated from the experimental data
    %   FSR_range = range that creates the acceptance tolerance.
    %   Inputted as a decimal (i.e., 0.5 = 50%)
    %   x = the concentration array for x
    %   t = time array from the model simulation
    %   eqDur = equilibrium duration
    %   simDur_end = t at the end of the simulation
    %   species = column in the x array that is being evalulated (e.g.,
    %   plasma leucine = 4)

    % FSR calculation
    mps_integral = cumtrapz(t(eqDur:simDur_end),...
        x(eqDur:simDur_end, species));
    mps_integral=mps_integral(end);
    LeucineMolarMass = 1/131.17; %mol/g
    FSR_grams_simulation = mps_integral/LeucineMolarMass; % MPS integral converted to FSR in grams of leucine

%     % Logical evaluation if the simulated value exists in the acceptance range.
%         % The range of the peak concentration is calculated using the 
%         % difference between the peak and basal values
%     rangeLow = FSR_grams_simulation*(1-FSR_range);
%     rangeHigh = FSR_grams_simulation*(1+FSR_range);
% 
%     pass = discretize(dataFSR, [rangeLow, rangeHigh])==1;

% UPDATED - above version does not use the simulated FSR value and dataFSR
% value correctly

    % Logical evaluation if the simulated value exists in the acceptance range.
        % The range of the peak concentration is calculated using the 
        % difference between the peak and basal values
    rangeLow = dataFSR*(1-FSR_range);
    rangeHigh = dataFSR*(1+FSR_range);

    pass = discretize(FSR_grams_simulation, [rangeLow, rangeHigh])==1;

end
