function [pass] = basalReturn_acceptanceCons(basalConc, peakConc, endConc_ratio, x, simDur_end, species)
    % Function assesses if the parameter set returns to basal levels at the
    % end of the MPSA simulation. *Conservative approach for acceptance*.
    % The acceptance is calculated as ratio (endConc_ratio) of the difference between the
    % peak and basal concentrations of the species
    %   endConc_ratio = ratio that creates the acceptance tolerance. I.e.,
    %   a factor of the difference between basal and peak concentration
    %   (e.g., 0 <-> (peak-basal)*ratio)
    %   x = the concentration array for x
    %   t = time array from the model simulation
    %   simDur_end = end of simulation (equilibrium+simulation)
    %   species = column in the x array that is being evalulated (e.g.,
    %   plasma leucine = 4)

    % concentration at the end of the simulation
    endConc_simulation = x(simDur_end, species);

    % Logical evaluation if the simulated value exists in the acceptance range.
    rangeHigh = (peakConc-basalConc)*endConc_ratio + basalConc;

    pass = discretize(endConc_simulation, [0, rangeHigh])==1;

end