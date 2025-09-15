function [pass] = peakConcentration_acceptance(peakConc, peakConc_range, basalConc, x, eqDur, species)
    % Function assesses the acceptance of the peak concentration of the 
    % MPSA simulation. 
    %   peakConc = Peak concentration of the data (or concentration around 
    %   which the range is calculated)
    %   peakConc_range = range that creates the acceptance tolerance. 
    %   Inputted as a decimal (i.e., 0.25 = 25%)
    %   basalConc = basal concentration of the data (or concentration
    %   around which the range is calcualted)
    %   x = the concentration array for x
    %   t = time array from the model simulation
    %   eqDur = equilibrium duration
    %   species = column in the x array that is being evalulated (e.g.,
    %   plasma leucine = 4)

    % simulated peak concentration following the end of the model
    % equilibration
    peakConc_simulation = max(x(eqDur:end, species));
    
    % Logical evaluation if the simulated value exists in the acceptance range.
        % The range of the peak concentration is calculated using the 
        % difference between the peak and basal values
    rangeLow = peakConc-(peakConc-basalConc)*(1-peakConc_range);
    rangeHigh = peakConc+(peakConc-basalConc)*(1-peakConc_range);
    
    pass = discretize(peakConc_simulation, [rangeLow, rangeHigh])==1;

end