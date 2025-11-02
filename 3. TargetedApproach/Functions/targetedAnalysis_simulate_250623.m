function [MPS_simulation, MPB_simulation, NB_simulation, mps_basal, mpb_basal, maxToBasal_Akt_pp, pToTotal_p70S6K_0min, pToTotal_p70S6K_30min, pToTotal_p70S6K_3hr] = ...
    targetedAnalysis_simulate_250623(subjectData, leucineInput, concentrationsIV, kValues, eqDur, simDur, ...
    expDuration, glucInfRate, leuInfRate, FilePathCumulative, p70S6K_basalActivity, insulinResistance)


    %% Model inputs
    AminoAcidInput = leucineInput; % grams
    LeucineMolarMass = 1/131.17; %mol/g
    AminoAcidInput_mass = AminoAcidInput*LeucineMolarMass;
    expDur = expDuration; %minutes  
    %% Calculating values for model simulation
    plasmaVolume = plasmaVolumeFxn(subjectData(1), subjectData(2), subjectData(3));
    skeletalMuscleVolume = skeletalMuscleVolumeFxn(subjectData(1), subjectData(2), subjectData(4), subjectData(5));

    %% Running model simulation
%     for k = 1:trials
%         count = k; disp(count)
        
        % selecting concentration parameter sets and conversion to moles
%         x0 = concentrationsIV(:,k);
        x0 = concentrationsIV;
        x0_mass = x0_massFxn(x0, plasmaVolume, skeletalMuscleVolume); 
        
%         kValues_k = kValues(:,k);
        kValues_k = kValues;

%         warning on verbose
%         s = warning('error', 'MATLAB:nearlySingularMatrix');

        try 
            % equilibrium period
            [t,x] = ode23s(@(t,x) modelDerivative_250623_targetedSim(t, x, kValues_k, [], ...
                [], [], eqDur, glucInfRate, leuInfRate, p70S6K_basalActivity, insulinResistance), ...
                0 :0.1: eqDur, x0_mass);

            % simulation (leucine feeding)
            x_length = size(x,1); %vertical length (rows)
            x0_2 = x(x_length, :);
            x0_2(1) = AminoAcidInput_mass;
    
            [t2,x2] = ode23s(@(t2,x2) modelDerivative_250623_targetedSim(t2, x2, kValues_k, [], ...
                [], [], eqDur, glucInfRate, leuInfRate, p70S6K_basalActivity, insulinResistance), ...
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

           %% MPS, MPB, NB Calculation
            eqDur_t = find(t==eqDur); eqDur_t = eqDur_t(1);  % locating the start of the simulation in the t vector
            expDur_end_t = find(t==(eqDur+expDur));  % locating the end of the simulation (180 minutes)
            LeucineMolarMass = 1/131.17; %mol/g

            %% FSR/MPS 
            MPS_species = 44;
            mps_integral = cumtrapz(t(eqDur_t:expDur_end_t), ...
                x(eqDur_t:expDur_end_t, MPS_species));
            mps_integral=mps_integral(end);
%             MPS_simulation(1,k) = mps_integral/LeucineMolarMass; % MPS integral converted to FSR in grams of leucine
            MPS_simulation = mps_integral/LeucineMolarMass; % MPS integral converted to FSR in grams of leucine
%             disp(MPS_simulation)

            %% MPB calculation
            MPB_species = 45; 
            mpb_integral = cumtrapz(t(eqDur_t:expDur_end_t), ...
                x(eqDur_t:expDur_end_t, MPB_species));
            mpb_integral = mpb_integral(end);
%             MPB_simulation(1,k) = mpb_integral/LeucineMolarMass; % MPS integral converted to FSR in grams of leucine
            MPB_simulation = mpb_integral/LeucineMolarMass; % MPS integral converted to FSR in grams of leucine
%             disp(MPB_simulation)

            %% NB calculation
%             NB_simulation(1,k) = (mps_integral/LeucineMolarMass) - (mpb_integral/LeucineMolarMass);
            NB_simulation = (mps_integral/LeucineMolarMass) - (mpb_integral/LeucineMolarMass);
%             disp(NB_simulation)


%% Anabolic resistance mechanisms
        
    % Insulin resistance
        % Peak p-Akt(S473,T308) during intervention
        maxAkt_pp = max(x_conc(3001:5402,22));
        basalAkt_pp = x_conc(3001,22);
        maxToBasal_Akt_pp = maxAkt_pp/basalAkt_pp;

    % mTORC1 insensitivity
        p70S6K_total = x_conc(1,27)+x_conc(1,28);
        p_p70S6K_30min = x_conc(3302,28);
        p_p70S6K_3hr = x_conc(4802,28);
        pToTotal_p70S6K_30min = p_p70S6K_30min/p70S6K_total;
        pToTotal_p70S6K_3hr = p_p70S6K_3hr/p70S6K_total;       
    
    % Elevated post-absorptive p-p70S6K
        p_p70S6K_t0 = x_conc(3001,28);
        total_p70S6K_t0 = x_conc(3001, 28)+x_conc(3001,27);
        pToTotal_p70S6K_0min = p_p70S6K_t0/total_p70S6K_t0;

%         % p70S6K signalling (plot)
%         plot(t, x_conc(:,28))
%         ylim([0, 9e-9])
%         grid on
%         hold on
%         plot(t, x_conc(:,27))
%         plot(t, (x_conc(:,27)+x_conc(:,28)))
%         hold off

%         % Insulin signalling modification (plot)
%         Akt_total = x_conc(:,19) + x_conc(:,20) + x_conc(:,21) + x_conc(:,22);
%         plot(t, Akt_total)
%         ylim([0, Akt_total(1,1)*1.1])
%         grid on
%         hold on
%         plot(t, x_conc(:,22))
%         legend({'Total Akt', 'p-Akt^{S473,T308}'});
%         hold off



        % Postabsorptive protein metabolism rates
        mps_basal = x(3001,44);
        mpb_basal = x(3001,45);

        catch
%             MPS_simulation(1,k) = nan;
            MPS_simulation = nan;
            disp('error')

            MPB_simulation = nan;
            NB_simulation = nan;
            p70S6K_oldFC = nan;
            mTORC1_update = nan; 
            mps_basal = nan; 
            mpb_basal = nan;
            maxToBasal_Akt_pp = nan;
            pToTotal_p70S6K_30min = nan;
            pToTotal_p70S6K_3hr = nan;
        end
end
