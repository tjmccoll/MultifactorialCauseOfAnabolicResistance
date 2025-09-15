function [pass] = peakConcentrationTiming_acceptance(peakTime, peakTime_range, x, t, eqDur, species)
    % Function assesses the acceptance of the peak concentration timing of
    % the MPSA simulation
    %   peakTime = Timing of peak concentration of the data (or 
    %   concentration around which the range is calculated)
    %   peakTime_range = range that creates the acceptance tolerance.
    %   Inputted as a two value list (i.e., [low, high])
    %   x = the concentration array for x
    %   t = time array from the model simulation
    %   eqDur = equilibrium duration
    %   species = column in the x array that is being evalulated (e.g.,
    %   plasma leucine = 4)

% simulated time of peak concentration following the end of the model
% equilibration
% *update* - added "sum(... ,2)" so that if two species are inputted, the
% columns are summed to allow for proper peak concentration and time 
peakConc_simulation = max(sum(x(eqDur:end, species),2));
peakConcTime_simulation = t(sum(x(eqDur:end, species),2) == peakConc_simulation);

% logical evaluation if the simulated value exists in the acceptance range
rangeLow = peakTime+peakTime_range(1);
rangeHigh = peakTime+peakTime_range(2);

pass = discretize(peakConcTime_simulation, [rangeLow, rangeHigh])==1;

end