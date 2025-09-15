function [pass] = basalConcentration_acceptance(basalConc, basalConc_range, x, eqDur, species)
    % Function assesses the acceptance of the basal concentration of the
    % MPSA simulation
    %   basalConc = basal concentration of the data (or concentration
    %   around which the range is calcualted)
    %   basalConc_range = range that creates the acceptance tolerance.
    %   Inputted as a decimal (i.e., 0.25 = 25%)
    %   x = the concentration array for x
    %   eqDur = equilibrium duration
    %   species = column in the x array that is being evalulated (e.g.,
    %   plasma leucine = 4)

% simulated basal concentration at the end of the model equilibrium
basal_simulation = x(eqDur, species);

% logical evaluation if the simulated value exists in the acceptance range
rangeLow = basalConc*(1-basalConc_range);
rangeHigh = basalConc*(1+basalConc_range);

pass = discretize(basal_simulation, [rangeLow, rangeHigh])==1;

end