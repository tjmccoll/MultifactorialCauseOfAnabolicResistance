function [passTest] = MPSArun_231006(trials, subjectData, concentrationsIV, kValues, eqDur, simDur, glucInfRate, leuInfRate, FilePathCumulative)
    % Function to run the model ODE, assess the acceptance of the model
    % simulation, and output the required MPSA data
    %   Trials = the number of MPSA trials
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

    %% pre-defining outputs
    passTest = nan(1,trials);

    %% Model inputs
    AminoAcidInput = 3.5; % grams
    LeucineMolarMass = 1/131.17; %mol/g
    AminoAcidInput_mass = AminoAcidInput*LeucineMolarMass;

    %% Calculating values for model simulation
    plasmaVolume = plasmaVolumeFxn(subjectData(1), subjectData(2), subjectData(3));
    skeletalMuscleVolume = skeletalMuscleVolumeFxn(subjectData(1), subjectData(2), subjectData(4), subjectData(5));

    %% Running model simulation

    for k = 1:trials
        count = k; disp(count)
        
        % selecting concentration parameter sets and conversion to moles
        x0 = concentrationsIV(:,k);
        x0_mass = x0_massFxn(x0, plasmaVolume, skeletalMuscleVolume); 
        
        kValues_k = kValues(:,k);

        warning on verbose
        s = warning('error', 'MATLAB:nearlySingularMatrix');

            try 
                % equilibrium period
                [t,x] = ode23s(@(t,x) modelDerivative_OIM_230922(t, x, kValues_k, [], ...
                    [], [], eqDur, glucInfRate, leuInfRate), ...
                    0 :0.1: eqDur, x0_mass);
    
                % simulation (leucine feeding)
                x_length = size(x,1); %vertical length (rows)
                x0_2 = x(x_length, :);
                x0_2(1) = AminoAcidInput_mass;
        
                [t2,x2] = ode23s(@(t2,x2) modelDerivative_OIM_230922(t2, x2, kValues_k, [], ...
                    [], [], eqDur, glucInfRate, leuInfRate), ...
                    eqDur :0.1: (eqDur+simDur), x0_2);
    
                x_bolus = [x; x2];
                t_bolus = [t; t2];
                
                    % Fm,a
                    x_bolus(:,41) = kValues_k(6)*x_bolus(:,4); % Fm,a
                    % Fm,0
                    x_bolus(:,42) = kValues_k(9)*x_bolus(:,8)./x_bolus(:,20) + ... % Akt(T)-mediated MPB
                        kValues_k(69)*x_bolus(:,8)./x_bolus(:,26) + ... % mTORC1-mediated MPB
                        kValues_k(11)*x_bolus(:,7) ; % KIC reamination
                    % Net Balance
                    x_bolus(:,43) = kValues_k(6)*x_bolus(:,4) ...
                        - kValues_k(7)*x_bolus(:,6);
                    % MPS
                    x_bolus(:,44) = kValues_k(15)*x_bolus(:,6).*x_bolus(:,28); % MPS, r15
                    % MPB
                    x_bolus(:,45) = kValues_k(9)*x_bolus(:,8) ./ x_bolus(:,20) + ... % Akt(T)-mediated MPB
                        kValues_k(69)*x_bolus(:,8) ./ x_bolus(:,26); % mTORC1-mediated MPB
                  
                    % parameter 46 set as nan's to allow for cost calculation to work
                    % (the x_bolus array needs to contain values for each cost parameter to be calculated)
                    x_bolus(:,46) = nan;
                    
                x = x_bolus;
                t = t_bolus;

                % Converting moles to concentration
                x_conc = x_concFxn(x, x0_mass, plasmaVolume, skeletalMuscleVolume);

                MatrixSingular=0;
    
                        % Assessing acceptibility criteria if parameter set
                        % runs
                        passTest(1,k) = MPSAcriteria_231006(t, x, x_conc, eqDur);
                        numPass =  sum(passTest==1);
                        disp(['Number of passed sets: ',num2str(numPass)]) % live output of the number of no passing parameter sets
    
                        % Single plot of the time-course simulations of
                        % passing paramter sets
                        if passTest(1,k) == 1
                            plot_mpsaPass_sp(1,t,x_conc, x, simDur, eqDur, kValues_k)
                        end

            catch
                t=nan;
                x=nan;
                passTest(1,k)=0;
                MatrixSingular=1;
                PosReal=0;

                disp('error')
    
                s;       %outputs error
                s('');  %resets error message
            end

    end

    % Save mpsa figure
    fileName_plot = sprintf('MPSA_passTest_plot - %s.eps', datestr(now, 'yy-mm-dd HH-MM-SS'));
    fileNameFull_plot = fullfile(FilePathCumulative, fileName_plot);
    
    figure1 = figure(1);
    saveas(figure1, fileNameFull_plot, 'epsc');


end