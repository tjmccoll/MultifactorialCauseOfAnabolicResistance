function [passTest] = MPSAcriteria_231006(t, x, x_conc,eqDur)
    % Function to assess the acceptability of each MPSA simulation
    %   t = time array from the ODE simulation
    %   x = x array in moles from the ODE simulation
    %   x_conc = x_concentration array from the ODE simulation
    %   eqDur = equilibrium duration used to attain steady state
    %   passTest = an array that reports if a simulation is accepted (1) or
    %   is rejected (0)

% locating the simulation start
eqDur_t = find(t==eqDur); eqDur_t = eqDur_t(1);  

% locating the end of the simulation (180 minutes)
simDur_end = 180; % minutes
simDur_end_t = find(t==(eqDur+simDur_end)); 

% location of species within the x (concentration) array
plasmaIns = 3; plasmaLeu = 4; intraLeu = 6; F_ma = 41; F_m0 = 42;
pAkt_S = [21,22]; % both species
p70S6K_T = 28; MPS = 44;


%% Qualitative assessment of plasma leucine
% Basal concentration, peak concentration and timing, return to basal
% levels

% Basal concentration
plLeu_basal = 1.23e-4; % mol/L; basal concentration from Glynn 2010
plLeu_basal_pass = basalConcentration_acceptance(plLeu_basal, 0.25, x_conc, eqDur_t, plasmaLeu);

% Peak concentration timing
plLeu_peakTime = 45; % minutes; Glynn 2010
plLeu_peakTime_pass = peakConcentrationTiming_acceptance(plLeu_peakTime, [-15,30], x_conc, t, eqDur_t, plasmaLeu);

% Peak concentration
plLeu_peakConc = 8.9e-4; % mol/L; Glynn 2010
plLeu_peakConc_pass = peakConcentration_acceptance(plLeu_peakConc, 0.5, plLeu_basal, x_conc, eqDur_t, plasmaLeu);

% Return to basal (assessment at end of simulation)
plLeu_endConc = 2.31e-4; % mol/L; Glynn 2010
plLeu_endConc_pass = basalReturn_acceptance(plLeu_endConc, 2, x_conc, simDur_end_t, plasmaLeu);

% Test that all species tests pass. If all tests=1, plLeu_pass=1 (i.e., passes)
plLeu_pass = all([plLeu_basal_pass, plLeu_peakTime_pass, plLeu_peakConc_pass, plLeu_endConc_pass]);


%% Qualitative assessment of plasma insulin
% Basal concentration, peak concentration and timing, return to basal
% levels

% Basal concentration
plIns_basal = 3.3e-11; % mol/L; basal concentration from Glynn 2010
plIns_basal_pass = basalConcentration_acceptance(plIns_basal, 0.25, x_conc, eqDur_t, plasmaIns);

% Peak concentration timing
plIns_peakTime = 30; % minutes; Glynn 2010
plIns_peakTime_pass = peakConcentrationTiming_acceptance(plIns_peakTime, [-15,30], x_conc, t, eqDur_t, plasmaIns);

% Peak concentration
plIns_peakConc = 9.2e-11; % mol/L; Glynn 2010
plIns_peakConc_pass = peakConcentration_acceptance(plIns_peakConc, 0.5, plIns_basal, x_conc, eqDur_t, plasmaIns);

% Return to basal (assessment at end of simulation)
plIns_endConc_pass = basalReturn_acceptanceCons(plIns_basal, plIns_peakConc, 0.5, x_conc, simDur_end_t, plasmaIns);

% Test that all species tests pass. If all tests=1, plLeu_pass=1 (i.e., passes)
plIns_pass = all([plIns_basal_pass, plIns_peakTime_pass, plIns_peakConc_pass, plIns_endConc_pass]);


%% Qualitative assessment of intracellular leucine
% Basal concentration, peak concentration and timing

% Basal concentration
intLeu_basal = 1.28e-4; % mol/L; basal concentration from Drummond 2010
intLeu_basal_pass = basalConcentration_acceptance(intLeu_basal, 0.5, x_conc, eqDur_t, intraLeu);

% Peak concentration timing
intLeu_peakTime = 120; % minutes; Drummond 2010
intLeu_peakTime_pass = peakConcentrationTiming_acceptance(intLeu_peakTime, [-60,30], x_conc, t, eqDur_t, intraLeu);

% Peak concentration
intLeu_peakConc = 2.8e-4; % mol/L; Drummond 2010
intLeu_peakConc_pass = peakConcentration_acceptance(intLeu_peakConc, 0.5, intLeu_basal, x_conc, eqDur_t, intraLeu);

% Test that all species tests pass. If all tests=1, plLeu_pass=1 (i.e., passes)
intLeu_pass = all([intLeu_basal_pass, intLeu_peakTime_pass, intLeu_peakConc_pass]);


%% Qualitative assessment of Fm,a
% Basal concentration, peak concentration and timing, return to basal
% levels

% Basal concentration
F_ma_basal = 2.29e-5; % mol/min; Glynn 2010
F_ma_basal_pass = basalConcentration_acceptance(F_ma_basal, 0.25, x, eqDur_t, F_ma);

% Peak concentration timing 
F_ma_peakTime = 45; % minutes, Glynn 2010
F_ma_peakTime_pass = peakConcentrationTiming_acceptance(F_ma_peakTime, [-15,30], x, t, eqDur_t, F_ma);

% Peak concentration
F_ma_peakConc = 1.28e-4; % mol/min; Glynn 2010
F_ma_peakConc_pass = peakConcentration_acceptance(F_ma_peakConc, 0.5, F_ma_basal, x, eqDur_t, F_ma);

% Return to basal (assessment at end of simulation)
F_ma_endConc_pass = basalReturn_acceptanceCons(F_ma_basal, F_ma_peakConc, 0.5, x, simDur_end_t, F_ma);

% Test that all species tests pass. If all tests=1, plLeu_pass=1 (i.e., passes)
F_ma_pass = all([F_ma_basal_pass, F_ma_peakTime_pass, F_ma_peakConc_pass, F_ma_endConc_pass]);


%% Qualitative assessment of Fm,0
% not currently included 


%% Qualitative assessment of p-Akt(S)
% basal phospho-to-total ratio, peak phospho-data and timing

% Basal phospho-to-total ratio
% Akt_totalSpecies = [19:22];
Akt_nonPhos_species = 19;

% pAkt_S_basal_acceptRange = [2, 50]; % percent
pAkt_S_basal_acceptRange = [0.5, 50]; % percent
% pAkt_S_basal_pass = basalPhosphoRatio_acceptance(pAkt_S_basal_acceptRange, x_conc, eqDur_t, pAkt_S, Akt_totalSpecies);
pAkt_S_basal_pass = basalPhosphoRatio_acceptance_231006(pAkt_S_basal_acceptRange, x_conc, eqDur_t, pAkt_S, Akt_nonPhos_species);

% Peak concentration timing
pAkt_S_peakTime = 60; % Scoping review
pAkt_S_peakTime_pass = peakConcentrationTiming_acceptance(pAkt_S_peakTime, [-15,60], x_conc, t, eqDur_t, pAkt_S);

% Peak concentration (fold change)
pAkt_S_peakFC_acceptRange = [1.1, 2.5]; %fold change range
pAkt_S_peakFC_pass = peakConcentrationFC_acceptance(pAkt_S_peakFC_acceptRange, x_conc, eqDur_t, pAkt_S);

% Test that all species tests pass. If all tests=1, plLeu_pass=1 (i.e., passes)
pAkt_S_pass = all([pAkt_S_basal_pass, pAkt_S_peakTime_pass, pAkt_S_peakFC_pass]);


%% Qualitative assessment of p-p70S6K(T)
% basal phospho-to-total ratio, peak phospho-data and timing

% Basal phospho-to-total ratio
% p70S6K_totalSpecies = [27:28];
p70S6K_nonPhos_species = 27;

p70S6K_T_basal_acceptRange = [1,40]; % percent
% p70S6K_T_basal_pass = basalPhosphoRatio_acceptance(p70S6K_T_basal_acceptRange, x_conc, eqDur_t, p70S6K_T, p70S6K_totalSpecies);
p70S6K_T_basal_pass = basalPhosphoRatio_acceptance_231006(p70S6K_T_basal_acceptRange, x_conc, eqDur_t, p70S6K_T, p70S6K_nonPhos_species);

% Peak concentration timing
p70S6K_T_peakTime = 90; % minutes, scoping review
p70S6K_T_peakTime_pass = peakConcentrationTiming_acceptance(p70S6K_T_peakTime, [-45, 60], x_conc, t, eqDur_t, p70S6K_T);

% Peak concentration (fold change)
p70S6K_T_peakFC_acceptRange = [1.25, 3]; % fold change range
p70S6K_T_peakFC_pass = peakConcentrationFC_acceptance(p70S6K_T_peakFC_acceptRange, x_conc, eqDur_t, p70S6K_T);

% Test that all species tests pass. If all tests=1, plLeu_pass=1 (i.e., passes)
p70S6K_T_pass = all([p70S6K_T_basal_pass, p70S6K_T_peakTime_pass, p70S6K_T_peakFC_pass]);


%% Qualitative assessment of FSR and MPS
% FSR integral, MPS peak timing

% FSR integral
FSR_integral = 0.36; % grams; Glynn 2010
FSR_integral_pass = FSR_acceptance(FSR_integral, 0.5, t, x, eqDur_t, simDur_end_t, MPS);

% Peak MPS timing 
MPS_peakTime = 60; % minutes, Glynn 2010
MPS_peakTime_pass = peakConcentrationTiming_acceptance(MPS_peakTime, [-15,60], x, t, eqDur_t, MPS);

% Test that all species tests pass. If all tests=1, plLeu_pass=1 (i.e., passes)
FSR_MPS_pass = all([FSR_integral_pass, MPS_peakTime_pass]);


%% Assessing if the criteria for all species passes

passTest = all([plLeu_pass, intLeu_pass, plIns_pass, F_ma_pass, pAkt_S_pass, p70S6K_T_pass, FSR_MPS_pass]);
passTest = double(passTest); %convert logical to numeric value

end



