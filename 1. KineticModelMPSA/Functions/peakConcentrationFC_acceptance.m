function [pass] = peakConcentrationFC_acceptance(peakFC_range, x, eqDur, species)
    % Function assesses the acceptance of the peak concentration of the 
    % MPSA simulation using fold change (phospho-data specific). 
    %   peakFC = Peak fold change of the data (or fold change around 
    %   which the range is calculated)
    %   peakFC_range = range that creates the acceptance tolerance. 
    %   Inputted as a two value list (i.e., [1,2] - 1 to 2 fold change)
    %   x = the concentration array for x
    %   t = time array from the model simulation
    %   eqDur = equilibrium duration
    %   species = column in the x array that is being evalulated (e.g.,
    %   plasma leucine = 4)

    % simulated basal concentration (end of steady state)
    basalConc_simulation = sum(x(eqDur, species));

    % simulated peak concentration following the end of the model
    % equilibration
    peakConc_simulation = max(sum(x(eqDur:end, species),2));
    
    % Peak concentration fold change
    peakFC_simulation = peakConc_simulation/basalConc_simulation;

    % Logical evaluation if the simulated FC exists in the acceptance FC range.
    pass = discretize(peakFC_simulation, [peakFC_range(1), peakFC_range(2)])==1;
    
end