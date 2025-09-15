function [pass] = basalPhosphoRatio_acceptance(phosphoRatioRange, x, eqDur, species, speciesTotal)
    % Function assesses the acceptance of the basal concentration of the
    % MPSA simulation
    %   phosphoRatioRange = rrange that creates the acceptance tolerance. 
    %   Inputted as a two value list (i.e., [1,40] - 1-40% phospho-to-total
    %   protein)
    %   x = the concentration array for x
    %   eqDur = equilibrium duration
    %   species = column in the x array that is being evalulated (e.g.,
    %   plasma leucine = 4)
    %   speciesTotal = columns in the x array that contain all
    %   phosphorylated and non-phosphorylated species of the protein

    % simulated basal concentration at the end of the model equilibrium
    basalPhospho_simulation = x(eqDur, species);

    % simulated total protein concentration at the end of model equilibrium
    basalTotal_simulation = sum(x(eqDur, speciesTotal));

    % phospho-to-total protein ratio
    phosphoToTotal_ratio = basalPhospho_simulation/basalTotal_simulation*100;
    
    % logical evaluation if the simulated value exists in the acceptance range    
    pass = discretize(phosphoToTotal_ratio, ...
        [phosphoRatioRange(1), phosphoRatioRange(2)])==1;

end