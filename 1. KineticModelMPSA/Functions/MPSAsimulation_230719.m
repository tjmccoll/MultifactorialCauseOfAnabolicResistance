function [t, x, x_conc] = MPSAsimulation_230719(subjectData, concentrationsIV, kValues, eqDur, simDur, glucInfRate, leuInfRate)
    % Function to run the model ODE with the accepted parameter sets
    %   SubjectData = subject specific data for calculating plasma and 
    %   skeletal muscle volume. SubjectData = [gender, height, mass, age, 
    %   bia]
    %   concentrationsIV = vector of initial concentrations for model 
    %   species
    %   kValues = inputted vector of rate constants
    %   eqDur = equilibrium duration to allow for model steady state (min)
    %   simDur = duration of model simulation (min)
    %   GlucoseInfRate = Glucose infusion rate required for the Sturis 
    %   module
    %   leuInfRate = Leucine infusion rate required to maintain leucine at 
    %   a steady state during the equilibrium period

    %% Model inputs
    AminoAcidInput = 3.5; % grams
    LeucineMolarMass = 1/131.17; %mol/g
    AminoAcidInput_mass = AminoAcidInput*LeucineMolarMass;

    %% Calculating values for model simulation
    plasmaVolume = plasmaVolumeFxn(subjectData(1), subjectData(2), subjectData(3));
    skeletalMuscleVolume = skeletalMuscleVolumeFxn(subjectData(1), subjectData(2), subjectData(4), subjectData(5));

    %% Running model simulation
    x0 = concentrationsIV;
    x0_mass = x0_massFxn(x0, plasmaVolume, skeletalMuscleVolume); 

    % Equilibrium
    [t,x] = ode23s(@(t,x) modelDerivative_OIM_220719(t, x, kValues, [], ...
        [], [], eqDur, glucInfRate, leuInfRate), ...
        0 :0.1: eqDur, x0_mass);

    % Simulation (leucine feeding)
    x_length = size(x,1); %vertical length (rows)
    x0_2 = x(x_length, :);
    x0_2(1) = AminoAcidInput_mass;

    [t2,x2] = ode23s(@(t2,x2) modelDerivative_OIM_220719(t2, x2, kValues, [], ...
        [], [], eqDur, glucInfRate, leuInfRate), ...
        eqDur :0.1: (eqDur+simDur), x0_2);

    x_bolus = [x; x2];
    t_bolus = [t; t2];

        x_bolus(:,41) = kValues(6)*x_bolus(:,4); % Fm,a
        x_bolus(:,42) = kValues(9)*x_bolus(:,8)./x_bolus(:,10) + kValues(11)*x_bolus(:,7); % Fm,0; insulin mediated MPB
        x_bolus(:,43) = kValues(6)*x_bolus(:,4) - kValues(7)*x_bolus(:,6); % Net Balance
        x_bolus(:,44) = kValues(15)*x_bolus(:,6).*x_bolus(:,28); % MPS, r15
        x_bolus(:,45) = kValues(9)*x_bolus(:,8) ./ x_bolus(:,10); % MPB, r9 (p-IR concentration)

        % parameter 46 set as nan's to allow for cost calculation to work
        % (the x_bolus array needs to contain values for each cost parameter to be calculated)
        x_bolus(:,46) = nan;

    x = x_bolus;
    t = t_bolus;

    % Converting moles to concentration
    x_conc = x_concFxn(x, x0_mass, plasmaVolume, skeletalMuscleVolume);


end