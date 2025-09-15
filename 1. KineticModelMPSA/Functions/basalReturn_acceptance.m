function [pass] = basalReturn_acceptance(endConc, endConc_range, x, simDur_end, species)
    % Function assesses if the paramter sets returns to basal levels at the
    % end of the MPSA simulation
    %   endConc = concentration at the end of model simulation (180
    %   minutes)
    %   endConc_range = range that creates the acceptance tolerance. 
    %   Inputted as a fold change (i.e., 2x)
    %   x = the concentration array for x
    %   t = time array from the model simulation
    %   eqDur = equilibrium duration
    %   species = column in the x array that is being evalulated (e.g.,
    %   plasma leucine = 4)

    % concentration at the end of the simulation
    endConc_simulation = x(simDur_end, species);

    % Logical evaluation if the simulated value exists in the acceptance range.
    rangeLow = endConc*(1/endConc_range);
    rangeHigh = endConc*endConc_range;

    pass = discretize(endConc_simulation, [rangeLow, rangeHigh])==1;

end